"""把共享入口验证的 counts h5ad 转为 scTour 输入；不训练或产生疾病结论。"""
import argparse
from pathlib import Path
import anndata as ad
import scanpy as sc
from prepare_geo import validate_counts

def prepare(input_path, output, celltype_key, batch_key, celltypes, hvg):
    source = ad.read_h5ad(input_path)
    validate_counts(source)
    assert celltype_key in source.obs and batch_key in source.obs, 'Provide actual annotation/batch column names'
    assert set(celltypes).issubset(set(source.obs[celltype_key])), 'Requested cell type absent'
    data = source[source.obs[celltype_key].isin(celltypes)].copy()
    validate_counts(data)
    assert not data.obs[celltype_key].isna().any() and not data.obs[batch_key].isna().any()
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
