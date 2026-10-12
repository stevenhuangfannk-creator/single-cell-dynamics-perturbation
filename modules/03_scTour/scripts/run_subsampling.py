"""论文 Figure 2h/i：真实训练子集、全部细胞推断、方向对齐后的排序比较。"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

def verify_old_cached_input(input_path, inferred_path):
    # Older runs predate input SHA recording. Compare IDs and full CSR counts without loading Scanpy/torch.
    import h5py
    import numpy as np
    with h5py.File(input_path, 'r') as source, h5py.File(inferred_path, 'r') as result:
        for axis in ['obs', 'var']:
            left = source[axis][source[axis].attrs['_index']][:]
            right = result[axis][result[axis].attrs['_index']][:]
            assert np.array_equal(left, right), f'Cached {axis} IDs/order differ'
        def obs_values(file, key):
            value = file['obs'][key]
            return value['categories'][:][value['codes'][:]] if isinstance(value, h5py.Group) else value[:]
        for key in ['celltype', 'batch', 'n_genes_by_counts']:
            assert np.array_equal(obs_values(source, key), obs_values(result, key)), f'Cached annotation/QC differs: {key}'
        for name in ['data', 'indices', 'indptr']:
            left, right = source['X'][name], result['X'][name]
            assert left.shape == right.shape, 'Cached counts differ'
            for start in range(0, len(left), 1000000):
                assert np.array_equal(left[start:start+1000000], right[start:start+1000000]), 'Cached counts differ'

def subsample(config, input_path, output, fractions):
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    expected = json.loads(Path(config).read_text(encoding='utf-8'))
    assert expected['nepoch'] == 400, 'Formal subsampling requires 400 epochs'
    assert 95 in fractions or (out/'B_95'/'run_summary.json').exists(), 'Include 95% or provide a completed B_95 reference'
    assert all(0 < frac <= 95 for frac in fractions), 'Training percentages must be in (0, 95]'
    with Path(input_path).open('rb') as handle:
        input_hash = hashlib.file_digest(handle, 'sha256').hexdigest()
    for frac in dict.fromkeys([*fractions, 95]):
        folder = out / f'B_{frac:g}'
        if (folder / 'run_summary.json').exists():
            cached = json.loads((folder / 'run_summary.json').read_text(encoding='utf-8'))
            used = json.loads((folder/'config_used.json').read_text(encoding='utf-8'))
            assert cached['epochs'] == 400 and used['percent'] == frac / 100
            assert all(used[key] == value for key, value in expected.items() if key != 'percent'), 'Cached config differs; use a new output directory'
            if 'input_sha256' in cached:
                assert cached['input_sha256'] == input_hash, 'Cached input differs; use a new output directory'
            else:
                verify_old_cached_input(input_path, folder/'inferred.h5ad')
            continue
        with open(out / f'B_{frac:g}.log', 'w', encoding='utf-8') as log:
            subprocess.run([sys.executable, str(Path(__file__).with_name('run_sctour.py')),
                            '--config', config, '--input-path', input_path, '--output', str(folder),
                            '--percent', str(frac / 100)], stdout=log, stderr=subprocess.STDOUT, check=True)
    # Heavy numerical imports only after training children exit; keep the launcher small on 16-GB PCs.
    import numpy as np
    import pandas as pd
    from scipy.stats import spearmanr
    columns = ['cell_id', 'ptime_before_manual', 'ptime']
    reference = pd.read_csv(out / 'B_95' / 'inference.tsv.gz', sep='\t', usecols=columns,
                            dtype={'cell_id': str}, keep_default_na=False).set_index('cell_id')
    rows = []
    for frac in fractions:
        data = pd.read_csv(out / f'B_{frac:g}' / 'inference.tsv.gz', sep='\t', usecols=columns,
                           dtype={'cell_id': str}, keep_default_na=False).set_index('cell_id')
        assert np.array_equal(data.index, reference.index), 'cell order mismatch'
        raw = spearmanr(reference.ptime_before_manual, data.ptime_before_manual)
        aligned = spearmanr(reference.ptime, data.ptime)
        cfg = json.loads((out / f'B_{frac:g}' / 'run_summary.json').read_text(encoding='utf-8'))
        rows.append({'percent': frac, 'seed': cfg['random_state'], 'epochs': cfg['epochs'],
                     'n_cells_inferred': len(data), 'training_seconds': cfg['training_seconds'],
                     'manual_reverse': cfg['manual_reverse'], 'rho_before_alignment': raw.statistic,
                     'rho_aligned': aligned.statistic, 'pvalue_aligned': aligned.pvalue})
    pd.DataFrame(rows).to_csv(out / 'robustness.tsv', sep='\t', index=False)
    print(pd.DataFrame(rows).to_string(index=False), flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--config', required=True)
    p.add_argument('--input-path', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--fractions', type=float, nargs='+', default=[1, 5, 10, 20, 40, 60, 80, 95])
    subsample(**vars(p.parse_args()))
