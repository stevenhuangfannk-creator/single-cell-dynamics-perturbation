"""Download raw 10x counts, transfer published labels by barcode, and audit targets."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import tarfile
import urllib.request
import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse

SOURCES = {
    "pbmc3k_filtered_gene_bc_matrices.tar.gz": (
        "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",
        "847d6ebd9a1ec9a768f2be7e40ca42cbfe75ebeb6d76a4c24167041699dc28b5"),
    "pbmc3k_processed.h5ad": (
        "https://raw.githubusercontent.com/chanzuckerberg/cellxgene/main/example-dataset/pbmc3k.h5ad",
        "0db367b991dd95809732b218539ede489bea99113807f62ebd7ccc970025fe38"),
}
# Pre-specified transcription regulators; inclusion requires >=5% detection.
CANDIDATES = ["SPI1", "CEBPB", "IRF8", "IRF1", "STAT1", "JUN", "FOS", "REL", "NFKB1", "CEBPD"]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=Path("work/knk/data"))
    ap.add_argument("--module", type=Path, default=Path("modules/01_scTenifoldKnk"))
    ap.add_argument("--n-genes", type=int, default=600)
    args = ap.parse_args()
    args.data.mkdir(parents=True, exist_ok=True)
    manifest = []
    for filename, (url, expected) in SOURCES.items():
        path = args.data / filename
        if not path.exists():
            urllib.request.urlretrieve(url, path)
        sha = hashlib.sha256(path.read_bytes()).hexdigest()
        if sha != expected:
            raise ValueError(f"Checksum mismatch: {filename}")
        manifest.append(dict(file=filename, url=url, bytes=path.stat().st_size, sha256=sha))
    mex = args.data / "filtered_gene_bc_matrices/hg19"
    if not mex.exists():
        with tarfile.open(args.data / "pbmc3k_filtered_gene_bc_matrices.tar.gz") as tf:
            tf.extractall(args.data, filter="data")
    counts = io.mmread(mex / "matrix.mtx").tocsr()
    genes_raw = pd.read_csv(mex / "genes.tsv", sep="\t", header=None, names=["ensembl", "symbol"])
    cells_raw = pd.read_csv(mex / "barcodes.tsv", sep="\t", header=None)[0].astype(str)
    assert counts.shape == (len(genes_raw), len(cells_raw))
    assert np.isfinite(counts.data).all() and (counts.data >= 0).all()
    assert np.equal(counts.data, np.floor(counts.data)).all() and cells_raw.is_unique
    # Sum duplicate symbols, preserving first-occurrence order and gene IDs in provenance.
    symbols = pd.Index(genes_raw.symbol.astype(str).unique())
    row_map = symbols.get_indexer(genes_raw.symbol.astype(str))
    agg = sparse.csr_matrix((np.ones(len(row_map)), (row_map, np.arange(len(row_map)))), shape=(len(symbols), len(row_map)))
    counts = (agg @ counts).tocsr()
    labels = ad.read_h5ad(args.data / "pbmc3k_processed.h5ad")
    selected_cells = labels.obs_names[labels.obs["louvain"].astype(str) == "CD14+ Monocytes"]
    cell_idx = pd.Index(cells_raw).get_indexer(selected_cells)
    assert (cell_idx >= 0).all() and selected_cells.is_unique
    raw_cd14 = counts[:, cell_idx].tocsr()
    lib = np.asarray(raw_cd14.sum(axis=0)).ravel()
    mito = np.asarray([str(g).upper().startswith("MT-") for g in symbols])
    mito_frac = np.asarray(raw_cd14[mito].sum(axis=0)).ravel() / lib
    # Published cell QC is rechecked from original counts, not the scaled h5ad.X.
    valid_cells = (lib >= 500) & (mito_frac < 0.1)
    raw_cd14 = raw_cd14[:, valid_cells].tocsr()
    selected_cells = selected_cells[valid_cells]
    lib = lib[valid_cells]
    detect = np.asarray((raw_cd14 > 0).mean(axis=1)).ravel()
    log_cpm = raw_cd14.toarray() / lib[None, :] * 1e6
    log_cpm = np.log1p(log_cpm)
    variance = log_cpm.var(axis=1, ddof=1)
    excluded = np.asarray([bool(re.match(r"^(RP[LS][0-9]|RPLP[0-9]|RPSA|MRP[LS][0-9]|MT-|MTND[0-9]|MTCO[0-9]|MTATP[0-9]|MTCYB|MTRNR2L)", str(g).upper())) for g in symbols])
    valid_genes = (detect >= 0.05) & (variance > 0) & ~excluded
    ranked = np.flatnonzero(valid_genes)[np.argsort(-variance[valid_genes], kind="stable")]
    qualified = [g for g in CANDIDATES if g in symbols and valid_genes[symbols.get_loc(g)]]
    assert "SPI1" in qualified, "Pre-specified SPI1 fails expression QC"
    forced = [symbols.get_loc(g) for g in qualified]
    chosen = np.asarray(forced + [i for i in ranked if i not in forced][:args.n_genes-len(forced)])
    X = raw_cd14[chosen].astype(np.float64)
    objects = args.module / "objects"; results = args.module / "results"
    objects.mkdir(parents=True, exist_ok=True); results.mkdir(parents=True, exist_ok=True)
    io.mmwrite(objects / "pbmc_cd14_counts.mtx", X, field="integer")
    pd.Series(symbols[chosen]).to_csv(objects / "pbmc_cd14_genes.tsv", sep="\t", header=False, index=False)
    pd.Series(selected_cells).to_csv(objects / "pbmc_cd14_cells.tsv", sep="\t", header=False, index=False)
    meta = labels.obs.loc[selected_cells].copy()
    meta["original_count_total"] = lib
    meta.to_csv(results / "pbmc_cd14_metadata.csv", index_label="cell")
    metrics = pd.DataFrame(dict(gene=symbols, detection_rate=detect, mean_counts=np.asarray(raw_cd14.mean(axis=1)).ravel(), logCPM_variance=variance, excluded_ribo_mito=excluded, selected=np.isin(np.arange(len(symbols)), chosen)))
    metrics.to_csv(results / "gene_selection_audit.csv", index=False)
    metrics.set_index("gene").reindex(CANDIDATES).to_csv(results / "candidate_expression_audit.csv")
    np.savetxt(objects / "candidate_targets.txt", np.asarray(qualified), fmt="%s")
    prepared = ad.AnnData(X=X.T.tocsr(), obs=meta, var=pd.DataFrame(index=symbols[chosen]))
    prepared.layers["counts"] = prepared.X.copy()
    prepared.uns["provenance"] = dict(dataset="10x PBMC3k", annotation="Scanpy legacy tutorial / cellxgene", counts_source="original 10x MEX", selected_cell_type="CD14+ Monocytes")
    prepared.write_h5ad(objects / "pbmc_cd14_counts.h5ad")
    summary = dict(status="RUN_OK", sources=manifest, original_shape=[32738,2700], original_unique_gene_count=len(symbols), duplicate_symbols_collapsed=len(genes_raw)-len(symbols), annotated_cells=labels.n_obs, annotated_cd14_cells=len(cell_idx), selected_cells=len(selected_cells), selected_genes=len(chosen), eligible_gene_count=int(valid_genes.sum()), pre_specified_targets=CANDIDATES, qualified_targets=qualified, rules="Published CD14 labels by exact barcode; library>=500, mito<0.1; detection>=0.05; ribo/mito excluded; top logCPM variance plus qualified pre-specified regulators; no response-based gene/target selection")
    (results / "pbmc_input_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))

if __name__ == "__main__":
    main()
