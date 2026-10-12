"""把共享入口验证的 counts h5ad 转为 scTour 输入；不训练或产生疾病结论。"""
import argparse
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.io import mmread
from prepare_geo import validate_counts

def prepare(input_path, output, celltype_key, batch_key, celltypes, hvg):
    path = Path(input_path)
    if path.is_dir():
        # Exact plain MEX schema emitted by shared/prepare_input.py --export-mex.
        genes = pd.read_csv(path/'genes.tsv', sep='\t', header=None, dtype=str, keep_default_na=False)[0].to_numpy()
        cells = pd.read_csv(path/'barcodes.tsv', sep='\t', header=None, dtype=str, keep_default_na=False)[0].to_numpy()
        metadata = pd.read_csv(path/'metadata.tsv', sep='\t', dtype=str, keep_default_na=False)
        metadata = metadata.set_index(metadata.columns[0])
        assert metadata.index.is_unique and set(cells).issubset(metadata.index), 'MEX metadata IDs differ'
        matrix = mmread(path/'matrix.mtx').tocsr().T.tocsr()
        assert matrix.shape == (len(cells), len(genes)), 'MEX axes differ from IDs'
        source = ad.AnnData(matrix, obs=metadata.reindex(cells), var=pd.DataFrame(index=genes))
    else:
        source = ad.read_h5ad(path)
    validate_counts(source)
    # Official encoder weights are float32; cast only when all count values remain exact.
    if source.X.dtype != np.float32:
        converted = source.X.astype(np.float32)
        before = source.X.data if sparse.issparse(source.X) else source.X
        after = converted.data if sparse.issparse(converted) else converted
        assert np.array_equal(before, after.astype(before.dtype)), 'Counts are not exactly representable as float32'
        source.X = converted
    assert celltype_key in source.obs and batch_key in source.obs, 'Provide actual annotation/batch column names'
    assert set(celltypes).issubset(set(source.obs[celltype_key])), 'Requested cell type absent'
    data = source[source.obs[celltype_key].isin(celltypes)].copy()
    validate_counts(data)
    assert not data.obs[celltype_key].isna().any() and not data.obs[batch_key].isna().any()
    assert data.obs[batch_key].astype(str).str.strip().ne('').all(), 'Missing batch labels'
    data.obs['celltype'] = data.obs[celltype_key].astype('category')
    data.obs['batch'] = data.obs[batch_key].astype(str)
    sc.pp.calculate_qc_metrics(data, percent_top=None, log1p=False, inplace=True)
    sc.pp.filter_genes(data, min_cells=20)
    assert data.n_vars >= hvg, 'Fewer retained genes than requested HVGs'
    data.layers['counts'] = data.X.copy()
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    data.write_h5ad(out/'qc_counts.h5ad', compression='gzip')
    sc.pp.highly_variable_genes(data, flavor='seurat_v3', n_top_genes=hvg)
    data.var.to_csv(out/'gene_qc_hvg.tsv', sep='\t')
    data[:, data.var.highly_variable].copy().write_h5ad(out/'train_counts.h5ad', compression='gzip')
    print(f'Prepared {data.n_obs} cells, {hvg} HVGs; no model was trained')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-path', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--celltype-key', default='cell_type')
    parser.add_argument('--batch-key', default='sample')
    parser.add_argument('--celltypes', nargs='+', required=True)
    parser.add_argument('--hvg', type=int, default=2000)
    prepare(**vars(parser.parse_args()))
