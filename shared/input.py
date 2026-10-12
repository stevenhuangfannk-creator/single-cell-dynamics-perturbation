"""Read raw counts without silently changing identifiers or count values."""
from pathlib import Path
import json

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread, mmwrite


def validate_counts(adata, min_cells=1, targets=()):
    if adata.n_obs < min_cells or adata.n_vars < 1:
        raise ValueError(f'Insufficient input: {adata.n_obs} cells, {adata.n_vars} genes')
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError('Cell and gene identifiers must be unique; resolve duplicates at source')
    values = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X).ravel()
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError('Counts contain non-finite or negative values')
    if not np.allclose(values, np.round(values), rtol=0, atol=1e-6):
        raise ValueError('Counts are not integer-valued; normalized/log data are not raw counts')
    missing = sorted(set(targets) - set(adata.var_names))
    if missing:
        raise ValueError(f'Targets missing from gene identifiers: {missing}')
    totals = np.asarray(adata.X.sum(axis=1)).ravel()
    if (totals == 0).any():
        raise ValueError('Empty cells found; document QC exclusions before analysis')
    detected = np.asarray((adata.X > 0).sum(axis=0)).ravel()
    return {'cells': adata.n_obs, 'genes': adata.n_vars, 'nnz': int(np.count_nonzero(values)),
            'count_sum': float(values.sum()), 'zero_genes': int((detected == 0).sum()),
            'library_size_min': float(totals.min()), 'library_size_median': float(np.median(totals)),
            'library_size_max': float(totals.max()),
            'targets': {g: {'detected_cells': int(detected[adata.var_names.get_loc(g)]),
                            'detection_fraction': float(detected[adata.var_names.get_loc(g)] / adata.n_obs)} for g in targets}}


def read_input(path, counts_layer='counts', metadata=None, gene_id='symbol'):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix == '.h5ad':
        source = ad.read_h5ad(path)
        if counts_layer == 'X':
            matrix = source.X
        elif counts_layer in source.layers:
            matrix = source.layers[counts_layer]
        else:
            raise ValueError(f'Explicit counts layer {counts_layer!r} missing; inspect source before selecting X')
        result = ad.AnnData(X=matrix.copy(), obs=source.obs.copy(), var=source.var.copy())
    elif path.is_dir():
        def locate(stem):
            for name in (stem, stem + '.gz'):
                candidate = path / name
                if candidate.exists():
                    return candidate
            raise FileNotFoundError(path / stem)
        matrix = mmread(locate('matrix.mtx')).tocsr().T.tocsr()
        features_path = next((p for p in (path / 'features.tsv', path / 'features.tsv.gz', path / 'genes.tsv', path / 'genes.tsv.gz') if p.exists()), None)
        if features_path is None:
            raise FileNotFoundError('features.tsv(.gz) or genes.tsv(.gz) missing')
        features = pd.read_csv(features_path, sep='\t', header=None, dtype=str, keep_default_na=False)
        barcodes = pd.read_csv(locate('barcodes.tsv'), sep='\t', header=None, dtype=str, keep_default_na=False)
        if matrix.shape != (len(barcodes), len(features)):
            raise ValueError('10x axes do not match barcodes and features')
        var = pd.DataFrame({'gene_id': features.iloc[:, 0].astype(str).values}, index=features.iloc[:, 1 if gene_id == 'symbol' and features.shape[1] > 1 else 0].astype(str).values)
        result = ad.AnnData(matrix, obs=pd.DataFrame(index=barcodes.iloc[:, 0].astype(str).values), var=var)
        if features.shape[1] >= 3:
            result = result[:, features.iloc[:, 2].eq('Gene Expression').values].copy()
    else:
        raise ValueError('Supported inputs: h5ad or 10x/Seurat-export MEX directory')
    if metadata:
        meta = pd.read_csv(metadata, sep='\t', converters={0: str})
        meta = meta.set_index(meta.columns[0])
        if not meta.index.is_unique or not result.obs_names.isin(meta.index).all():
            raise ValueError('Metadata identifiers are duplicated or missing input cells')
        overlap = set(meta.columns) & set(result.obs.columns)
        if overlap:
            raise ValueError(f'Metadata would overwrite existing annotations: {sorted(overlap)}')
        result.obs = result.obs.join(meta.loc[result.obs_names])
    result.layers['counts'] = result.X.copy()
    return result


def assess_pi(adata, cell_type_key, condition_key, batch_key, donor_key, targets=()):
    report = validate_counts(adata, targets=targets)
    keys = [cell_type_key, condition_key, batch_key, donor_key]
    missing = [key for key in keys if key not in adata.obs]
    if missing:
        raise ValueError(f'Required annotation columns missing: {missing}')
    if adata.obs[keys].isna().any().any():
        raise ValueError('Required annotations contain missing values')
    report['cell_type_condition_counts'] = adata.obs.groupby([cell_type_key, condition_key], observed=True).size().reset_index(name='cells').to_dict('records')
    report['donors_by_condition'] = adata.obs.groupby(condition_key, observed=True)[donor_key].nunique().to_dict()
    report['batches_by_condition'] = pd.crosstab(adata.obs[condition_key], adata.obs[batch_key]).to_dict()
    report['trajectory_evidence'] = 'NOT_EVALUATED: PI/Healthy labels alone do not establish developmental time'
    report['scope'] = 'Input assessment only; no PI/Healthy biological conclusions'
    return report


def save_json(report, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')


def export_mex(adata, path):
    """Export audited counts as genes × cells with one-column IDs for R."""
    validate_counts(adata)
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    matrix = sparse.coo_matrix(adata.X.T)
    matrix.data = np.rint(matrix.data).astype(np.int64)
    mmwrite(path / 'matrix.mtx', matrix, field='integer')
    (path / 'genes.tsv').write_text('\n'.join(adata.var_names) + '\n', encoding='utf-8')
    (path / 'barcodes.tsv').write_text('\n'.join(adata.obs_names) + '\n', encoding='utf-8')
    adata.obs.to_csv(path / 'metadata.tsv', sep='\t')
