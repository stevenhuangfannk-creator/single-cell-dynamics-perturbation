"""实际 scTour 1.0.0 训练、检查点、推断和嵌入。模型训练与绘图独立。"""
import argparse
import importlib.metadata
import hashlib
import json
import os
from pathlib import Path
import platform
import time
import psutil

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr
from scipy import sparse
import sctour as sct
import torch
from prepare_geo import validate_counts

def run(config, input_path, output, smoke=False, percent=None, seed=None, load=None, cpu=False):
    cfg = json.loads(Path(config).read_text(encoding='utf-8'))
    assert cfg['loss_mode'] in ['nb', 'zinb'], 'This validated runner supports raw-count NB/ZINB only'
    if percent is not None:
        cfg['percent'] = percent
    if seed is not None:
        cfg['random_state'] = seed
    assert 0 < cfg['percent'] < 1 and cfg['nepoch'] > 0, 'Require a non-empty training/validation split and positive epochs'
    cfg['use_gpu_requested'] = not cpu
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    torch.set_num_threads(4)
    np.random.seed(cfg['random_state'])
    torch.manual_seed(cfg['random_state'])
    adata = ad.read_h5ad(input_path)
    if cfg['loss_mode'] in ['nb', 'zinb']:
        validate_counts(adata)
    required_obs = ['celltype', 'batch', 'n_genes_by_counts']
    missing = [key for key in required_obs if key not in adata.obs]
    assert not missing, f'Missing observation columns: {missing}'
    assert not adata.obs[required_obs].isna().any().any(), 'Missing annotations or QC values'
    adata.obs['celltype'] = adata.obs.celltype.astype('category')
    stage_order = cfg.get('stage_order', [])
    assert len(stage_order) == len(set(stage_order)), 'stage_order contains duplicates'
    assert adata.n_vars == cfg['hvg'], 'Input HVG count differs from config'
    assert 0 < cfg['pca_components'] < min(adata.shape), 'PCA components exceed the input dimensions'
    with Path(input_path).open('rb') as handle:
        input_hash = hashlib.file_digest(handle, 'sha256').hexdigest()
    if smoke:
        # Stratified subset preserves every lineage state; explicitly a software check.
        indices = np.concatenate([np.flatnonzero(adata.obs.celltype == c)[:40]
                                  for c in adata.obs.celltype.cat.categories])
        adata = adata[indices].copy()
        cfg['nepoch'] = 3
        cfg['percent'] = .8
        cfg['batch_size'] = 128
    resource = {'python': platform.python_version(), 'cpu_count': os.cpu_count(),
                'host_ram_bytes': psutil.virtual_memory().total,
                'available_ram_bytes_at_start': psutil.virtual_memory().available,
                'cuda_available': torch.cuda.is_available(), 'torch_cuda': torch.version.cuda,
                'gpu': torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
                'gpu_memory_bytes': torch.cuda.get_device_properties(0).total_memory if torch.cuda.is_available() else 0,
                'versions': {p: importlib.metadata.version(p) for p in
                             ['sctour', 'torch', 'scanpy', 'anndata', 'numpy', 'scipy', 'pandas', 'matplotlib', 'torchdiffeq']}}
    (out / 'config_used.json').write_text(json.dumps(cfg, indent=2), encoding='utf-8')
    (out / 'environment.json').write_text(json.dumps(resource, indent=2), encoding='utf-8')
    started = time.perf_counter()
    model_args = {k: cfg[k] for k in ['percent', 'nepoch', 'batch_size', 'random_state', 'loss_mode', 'n_latent', 'alpha_recon_lec', 'alpha_recon_lode']}
    model_args['use_gpu'] = not cpu
    trainer = sct.predict.load_model(load) if load else sct.train.Trainer(adata, **model_args)
    if load:
        saved = trainer.adata
        assert saved.obs_names.equals(adata.obs_names) and saved.var_names.equals(adata.var_names), 'Checkpoint input IDs/order differ'
        if sparse.issparse(saved.X) and sparse.issparse(adata.X):
            same_counts = (saved.X != adata.X).nnz == 0
        else:
            x_saved = saved.X.toarray() if sparse.issparse(saved.X) else saved.X
            x_input = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
            same_counts = np.array_equal(x_saved, x_input)
        assert same_counts, 'Checkpoint counts differ; --load is for the original input only'
        assert np.array_equal(saved.obs.n_genes_by_counts, adata.obs.n_genes_by_counts), 'Checkpoint QC values differ'
        for key in ['percent', 'nepoch', 'batch_size', 'random_state', 'loss_mode']:
            assert getattr(trainer, key) == cfg[key], f'Checkpoint parameter differs: {key}'
    if not load:
        trainer.train()
    elapsed = time.perf_counter() - started
    original_time = trainer.get_time()
    adata.obs['ptime_before_manual'] = original_time
    # Biological check uses published stage labels; reversal is applied only when the medians are reversed.
    medians = pd.Series(original_time).groupby(adata.obs.celltype.to_numpy(), observed=True).median()
    present = [s for s in stage_order if s in medians.index]
    direction_rho = float(spearmanr(np.arange(len(present)), medians.reindex(present)).statistic) if len(present) >= 2 else None
    before = medians.reindex(present).to_numpy()
    # A partially scrambled trajectory cannot be repaired by reversing all times.
    reverse = bool(direction_rho is not None and direction_rho < 0 and np.all(np.diff(before) < 0))
    adata.obs['ptime'] = sct.train.reverse_time(original_time) if reverse else original_time
    # Official latent API and official high-dimensional differential equation output.
    mix, z, zt = trainer.get_latentsp(alpha_z=cfg['alpha_z'], alpha_predz=cfg['alpha_predz'])
    adata.obsm['X_TNODE'] = mix
    adata.obsm['X_z'] = z
    adata.obsm['X_zt'] = zt
    adata.obsm['X_VF'] = trainer.get_vector_field(original_time, mix)
    # X_VF remains the unmodified official output; plot_vector_field(reverse=True) corrects display.
    adata.uns['manual_reverse'] = bool(reverse)
    adata.uns['automatic_time_reverse'] = bool(trainer.time_reverse)
    for name in ['X_TNODE', 'X_z', 'X_zt', 'X_VF']:
        assert adata.obsm[name].shape == (adata.n_obs, cfg['n_latent'])
        assert np.isfinite(adata.obsm[name]).all(), name
    assert np.isfinite(adata.obs.ptime).all() and adata.obs.ptime.between(0, 1).all()
    if not load:
        trainer.save_model(str(out), 'model')
        # Save exact sampled training/validation cells, not merely the intended percentage.
        # scTour MakeDataset does not retain IDs, so replay and verify actual selected loader tensors.
        # The split is deterministic; replay seed before two-round random sampling.
        from sctour.data import split_index
        np.random.seed(cfg['random_state'])
        train_idx, val_idx = split_index(adata.n_obs, cfg['percent'], trainer.val_frac)
        expected_train = np.log1p(adata.X[train_idx])
        expected_val = np.log1p(adata.X[val_idx])
        if hasattr(expected_train, 'toarray'):
            expected_train = expected_train.toarray()
            expected_val = expected_val.toarray()
        assert np.array_equal(trainer.train_dataset.data.numpy(), expected_train), 'replayed train IDs do not match actual tensors'
        assert np.array_equal(trainer.val_dataset.data.numpy(), expected_val), 'replayed validation IDs do not match actual tensors'
        del expected_train, expected_val
        pd.DataFrame({'cell_id': adata.obs_names[train_idx]}).to_csv(out / 'train_cells.tsv', sep='\t', index=False)
        pd.DataFrame({'cell_id': adata.obs_names[val_idx]}).to_csv(out / 'validation_cells.tsv', sep='\t', index=False)
        pd.DataFrame(trainer.log).to_csv(out / 'loss.tsv', sep='\t', index_label='epoch_zero_based')
        # Official save/load API roundtrip: deterministic time encoder output is identical.
        reloaded = sct.predict.load_model(str(out / 'model.pth'))
        delta = float(np.max(np.abs(reloaded.get_time() - original_time)))
        assert delta < 1e-6, f'reloaded ptime differs: {delta}'
    else:
        delta = None
    # Distinct counts and normalized plotting AnnData preserve raw NB input.
    plot = adata.copy()
    sc.pp.normalize_total(plot, target_sum=1e4)
    sc.pp.log1p(plot)
    sc.pp.scale(plot, max_value=10)
    sc.tl.pca(plot, n_comps=cfg['pca_components'], svd_solver='arpack', random_state=cfg['random_state'])
    sc.pp.neighbors(plot, n_pcs=cfg['pca_components'], n_neighbors=15, key_added='pca_neighbors', random_state=cfg['random_state'])
    sc.tl.umap(plot, neighbors_key='pca_neighbors', min_dist=.1, random_state=cfg['random_state'])
    adata.obsm['X_pca'] = plot.obsm['X_pca']
    adata.obsm['X_umap_pca'] = plot.obsm['X_umap'].copy()
    sc.pp.neighbors(adata, use_rep='X_TNODE', n_neighbors=15, key_added='latent_neighbors', random_state=cfg['random_state'])
    sc.tl.umap(adata, neighbors_key='latent_neighbors', min_dist=.1, random_state=cfg['random_state'])
    adata.obsm['X_umap_latent'] = adata.obsm['X_umap'].copy()
    adata.write_h5ad(out / 'inferred.h5ad', compression='gzip')
    coordinates = adata.obs[['celltype', 'batch', 'ptime_before_manual', 'ptime']].copy()
    for rep in ['X_umap_pca', 'X_umap_latent', 'X_TNODE', 'X_VF']:
        for i in range(adata.obsm[rep].shape[1]):
            coordinates[f'{rep}_{i+1}'] = adata.obsm[rep][:, i]
    coordinates.to_csv(out / 'inference.tsv.gz', sep='\t', index_label='cell_id')
    memory = psutil.Process().memory_info()
    resource.update({'actual_model_device': str(trainer.device),
                     'compute_scope': 'official encoder/decoder on model device; official ODE on CPU',
                     'peak_process_working_set_bytes': getattr(memory, 'peak_wset', None),
                     'process_rss_bytes_at_end': memory.rss,
                     'torch_peak_allocated_bytes': torch.cuda.max_memory_allocated() if trainer.device.type == 'cuda' else 0,
                     'total_run_seconds': time.perf_counter() - started})
    (out / 'environment.json').write_text(json.dumps(resource, indent=2), encoding='utf-8')
    after_medians = adata.obs.groupby('celltype', observed=True).ptime.median().reindex(present)
    stages_monotonic = bool(np.all(np.diff(after_medians) > 0)) if len(present) >= 2 else None
    summary = {'status': 'RUN_OK' if smoke else 'PARTIAL', 'smoke': smoke,
               'input_sha256': input_hash,
               'n_cells_inferred': adata.n_obs, 'hvg': adata.n_vars, 'percent': cfg['percent'],
               'epochs': trainer.nepoch, 'training_seconds': elapsed, 'random_state': cfg['random_state'],
               'manual_reverse': bool(reverse), 'automatic_time_reverse': bool(trainer.time_reverse),
               'developmental_median_rho_before': direction_rho,
               'developmental_direction_assessment': 'ASSESSED' if len(present) >= 2 else 'NOT_ASSESSED',
               'developmental_medians_strictly_increasing': stages_monotonic,
               'ptime_median_by_celltype': adata.obs.groupby('celltype', observed=True).ptime.median().to_dict(),
               'roundtrip_max_ptime_difference': delta,
               'final_train_loss': trainer.log['train_loss'][-1] if not load else None,
               'final_validation_loss': trainer.log['validation_loss'][-1] if not load else None,
               'reason': '完成真实计算；与原图逐项比较及稳健性评估决定最终复现状态。'}
    (out / 'run_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(summary, ensure_ascii=True, indent=2), flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--config', required=True)
    p.add_argument('--input-path', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--smoke', action='store_true')
    p.add_argument('--percent', type=float)
    p.add_argument('--seed', type=int)
    p.add_argument('--load')
    p.add_argument('--cpu', action='store_true')
    run(**vars(p.parse_args()))
