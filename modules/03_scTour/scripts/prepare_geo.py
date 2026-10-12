"""从 GEO 原始表恢复 counts 和原始细胞元数据，再严格筛选论文谱系。"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

SOURCES = {
    'A': [('GSE95315_10X_expression_data_v2.tab.gz', 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE95nnn/GSE95315/suppl/')],
    'B': [('GSE104323_10X_expression_data_V2.tab.gz', 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104323/suppl/'),
          ('GSE104323_metadata_barcodes_24185cells.txt.gz', 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104323/suppl/')],
}
LABELS = {
    'A': ['nIPC', 'Neuroblast_1', 'Neuroblast_2', 'Granule-immature', 'Granule-mature'],
    'B': ['nIPC', 'Neuroblast', 'Immature-Pyr', 'Immature-GC', 'GC-juv', 'GC-adult'],
}
EXPECTED_SHA256 = {
    'GSE95315_10X_expression_data_v2.tab.gz': 'd2d2975854e07f1bf6ead0d57e94524de1ba987e83b165981b435c7901ec7aa6',
    'GSE104323_10X_expression_data_V2.tab.gz': '06046d4591e36a4b4bf4ceb6e9dd3212e639291e49438ec2d2500602241e9910',
    'GSE104323_metadata_barcodes_24185cells.txt.gz': '3e2587ceb8c556632462481328c0337fd3eaf2b683972de107818e60a8554ff4',
}

def validate_counts(adata):
    x = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X)
    assert adata.obs_names.is_unique and adata.var_names.is_unique, 'IDs must be unique'
    assert np.isfinite(x).all() and (x >= 0).all(), 'counts contain invalid values'
    assert np.equal(x, np.floor(x)).all(), 'NB requires integer counts'
    assert adata.n_obs > 20 and adata.n_vars > 20

def prepare(dataset, raw_dir, out_dir, hvg=None, download=False):
    raw_dir, out_dir = Path(raw_dir), Path(out_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for name, base in SOURCES[dataset]:
        p = raw_dir / name
        existing_hash = hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None
        if download and existing_hash != EXPECTED_SHA256[name]:
            subprocess.run(['curl', '-fL', '--retry', '5', '--retry-all-errors', '-C', '-', '-o', str(p), base + name], check=True)
        if not p.exists():
            raise FileNotFoundError(p)
        # Read through gzip EOF to verify CRC; SHA256 of the downloaded compressed file.
        with gzip.open(p, 'rb') as f:
            while f.read(1024 * 1024):
                pass
        actual_hash = hashlib.sha256(p.read_bytes()).hexdigest()
        assert actual_hash == EXPECTED_SHA256[name], f'Source changed or download damaged: {name}'
        records.append({'file': name, 'url': base + name, 'bytes': p.stat().st_size,
                        'sha256': actual_hash, 'gzip_crc': 'PASS'})
    expression = raw_dir / SOURCES[dataset][0][0]
    with gzip.open(expression, 'rt') as f:
        cells = next(f).rstrip('\r\n').split('\t')[1:]
        if dataset == 'A':
            ages = next(f).rstrip('\r\n').split('\t')[1:]
            labels = next(f).rstrip('\r\n').split('\t')[1:]
            meta = pd.DataFrame({'celltype_original': labels, 'age': ages}, index=cells)
            skiprows = [1, 2]
        else:
            meta = pd.read_csv(raw_dir / SOURCES['B'][1][0], sep='\t', index_col=0).dropna(how='all')
            meta = meta.reindex(cells)
            assert meta['characteristics: cell cluster'].notna().all(), 'metadata mismatch'
            meta = meta.rename(columns={'characteristics: cell cluster': 'celltype_original', 'characteristics: age': 'age'})
            skiprows = None
    meta['celltype'] = meta.celltype_original.replace({'nIPC-perin': 'nIPC'})
    meta['batch'] = ['_'.join(c.split('_')[:2]) for c in cells]
    select = meta.celltype.isin(LABELS[dataset]).to_numpy()
    meta = meta.loc[select].copy()
    meta['celltype'] = pd.Categorical(meta.celltype, categories=LABELS[dataset], ordered=True)
    # Read selected cell columns in bounded row chunks; retain the source order.
    usecols = [0] + (np.flatnonzero(select) + 1).tolist()
    matrices, genes = [], []
    reader = pd.read_csv(expression, sep='\t', skiprows=skiprows, usecols=usecols,
                         index_col=0, chunksize=256,
                         dtype={cells[i]: np.int32 for i in np.flatnonzero(select)})
    for chunk in reader:
        genes.extend(chunk.index.astype(str))
        matrices.append(sparse.csr_matrix(chunk.to_numpy(dtype=np.float32)))
    counts = sparse.vstack(matrices, format='csr').T.tocsr()
    adata = ad.AnnData(X=counts, obs=meta, var=pd.DataFrame(index=genes))
    validate_counts(adata)
    assert adata.n_obs == {'A': 4007, 'B': 15174}[dataset], 'paper lineage size mismatch'
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    all_genes = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=20)
    adata.layers['counts'] = adata.X.copy()
    adata.write_h5ad(out_dir / f'{dataset}_qc_counts.h5ad', compression='gzip')
    hvg = hvg or {'A': 500, 'B': 2000}[dataset]
    sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=hvg)
    adata.var.to_csv(out_dir / f'{dataset}_gene_qc_hvg.tsv', sep='\t')
    trained = adata[:, adata.var.highly_variable].copy()
    trained.write_h5ad(out_dir / f'{dataset}_train_counts.h5ad', compression='gzip')
    trained.obs.to_csv(out_dir / f'{dataset}_metadata.tsv', sep='\t')
    log = {'dataset': dataset, 'original_cells': len(cells), 'lineage_cells': adata.n_obs,
           'genes_before': all_genes, 'genes_after_min20': adata.n_vars, 'hvg': trained.n_vars,
           'celltype_counts': meta.celltype.value_counts(sort=False).to_dict(),
           'counts_validation': 'PASS', 'hvg_input': 'raw counts, seurat_v3', 'sources': records}
    (out_dir / f'{dataset}_preprocessing.json').write_text(json.dumps(log, indent=2), encoding='utf-8')
    print(json.dumps(log, indent=2), flush=True)
    return log

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--dataset', choices=['A', 'B'], required=True)
    p.add_argument('--raw-dir', required=True)
    p.add_argument('--out-dir', required=True)
    p.add_argument('--hvg', type=int)
    p.add_argument('--download', action='store_true')
    prepare(**vars(p.parse_args()))
