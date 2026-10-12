"""Run upstream README Python API on a supplied, audited real counts matrix."""
import argparse
import importlib.metadata
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.io import mmread
from virt_overexpression import scTenifoldOECUDA
import virt_knock.core as vk_core


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--genes", type=Path, required=True)
    parser.add_argument("--cells", type=Path, required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    torch.set_num_threads(2)
    vk_core.DEVICE = torch.device("cpu")
    gene_names = pd.read_csv(args.genes, sep="\t", header=None).iloc[:, 0].astype(str)
    cell_names = pd.read_csv(args.cells, sep="\t", header=None).iloc[:, 0].astype(str)
    counts = mmread(args.counts).toarray()
    if counts.shape != (len(gene_names), len(cell_names)):
        raise ValueError("Counts, gene IDs and cell IDs have different dimensions")
    if not np.isfinite(counts).all() or (counts < 0).any() or not np.equal(counts, np.floor(counts)).all():
        raise ValueError("Input must contain finite nonnegative integer counts")
    data = pd.DataFrame(counts, index=gene_names, columns=cell_names)
    started = time.perf_counter()
    # Published README API, with small bounded demonstration settings.
    model = scTenifoldOECUDA(data=data, oe_genes=[args.target], oe_fold=2,
                            n_nets=3, n_samp_cells=min(300, len(cell_names)),
                            n_comp=3, q=.9, ma_dim=2, random_state=1,
                            qc_kws={"min_lib_size": 0, "remove_outlier_cells": False,
                                    "min_percent": 0, "min_exp_sum": 0})
    result = model.run()
    args.output.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output / "native_README_api_dr.csv", index=False)
    model.wt_tensor.to_csv(args.output / "native_README_api_WT.csv")
    model.oe_tensor.to_csv(args.output / "native_README_api_OE.csv")
    row = model.wt_tensor.loc[args.target].to_numpy()
    assert (row > 0).any(), "Native demo entered fallback; inspect before interpretation"
    assert np.isfinite(result[["Distance", "FC", "adjusted_p-value"]]).all().all()
    metadata = {"status": "RUN_OK", "purpose": "Published README API execution; upstream provides no bundled reference dataset or numeric expected results",
                "data_type": "real counts; already selected/QC checked by Knk module",
                "genes": len(data), "cells": data.shape[1], "target": args.target,
                "n_nets": 3, "n_samp_cells": min(300, len(cell_names)), "ma_dim": 2,
                "oe_fold": 2, "seed": 1, "device": "cpu",
                "positive_target_row_edges": int((row > 0).sum()),
                "elapsed_seconds": time.perf_counter() - started,
                "versions": {p: importlib.metadata.version(p) for p in ["torch", "virt_knock", "virt_overexpression"]}}
    (args.output / "native_README_api_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
