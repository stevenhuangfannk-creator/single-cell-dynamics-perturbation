"""Experimental OE: audited upstream row amplification, official R comparisons."""
import argparse
import hashlib
import importlib.metadata
import json
import os
import subprocess
import time
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
from pathlib import Path

import numpy as np
import pandas as pd
MODULE = Path(__file__).resolve().parents[1]


def validate_network(wt, target):
    if wt.shape[0] != wt.shape[1] or not wt.index.equals(wt.columns):
        raise ValueError("WT must be square, with identical ordered row and column IDs")
    if not wt.index.is_unique or target not in wt.index:
        raise ValueError("Gene IDs must be unique and target must exist")
    if not np.isfinite(wt.to_numpy()).all() or not np.all(np.diag(wt) == 0):
        raise ValueError("WT must have finite values and zero diagonal")
    if not (wt.loc[target] > 0).any():
        raise ValueError("Target has no positive outgoing edges: upstream fallback is not accepted for this analysis")


def apply_oe(wt, target, fold):
    if int(pd.__version__.split(".")[0]) >= 3:
        raise RuntimeError("Upstream virt_overexpression requires pandas<3: homogeneous .values is read-only in pandas 3")
    if not np.isfinite(fold) or fold < 1:
        raise ValueError("OE fold must be finite and >=1")
    validate_network(wt, target)
    from virt_overexpression import virtual_overexpression
    oe = virtual_overexpression(wt, [target], oe_fold=fold, expression_data=None)
    expected = wt.copy()
    mask = expected.loc[target] > 0
    expected.loc[target, mask] *= fold
    np.testing.assert_array_equal(oe.to_numpy(), expected.to_numpy())
    return oe


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--wt", type=Path, required=True)
    parser.add_argument("--ko", type=Path, required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--folds", type=float, nargs="+", default=[1, 1.5, 2, 3])
    parser.add_argument("--rscript", required=True)
    parser.add_argument("--r-library", required=True)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--ma-dim", type=int, default=2)
    parser.add_argument("--max-genes", type=int, default=1000, help="Explicit resource limit; raise only after memory assessment")
    parser.add_argument("--output", type=Path, default=MODULE)
    args = parser.parse_args()
    wt = pd.read_csv(args.wt, index_col=0)
    ko = pd.read_csv(args.ko, index_col=0)
    validate_network(wt, args.target)
    if len(wt) > args.max_genes:
        raise ValueError(f"Network has {len(wt)} genes, exceeds --max-genes={args.max_genes}; assess O(genes²) memory first")
    if not wt.index.equals(ko.index) or not wt.columns.equals(ko.columns):
        raise ValueError("KO must preserve the exact WT ordered row and column IDs")
    expected_ko = wt.copy()
    expected_ko.loc[args.target] = 0
    np.testing.assert_array_equal(ko.to_numpy(), expected_ko.to_numpy())
    for folder in ["objects", "results", "reports"]:
        (args.output / folder).mkdir(parents=True, exist_ok=True)
    wt_path = args.output / "objects" / "WT.csv"
    wt.to_csv(wt_path)
    networks = [("KO", ko)]
    for fold in args.folds:
        networks.append((f"OE_{fold:g}", apply_oe(wt, args.target, fold)))
    metrics = []
    started = time.perf_counter()
    for name, network in networks:
        network_path = args.output / "objects" / f"{name}.csv"
        network.to_csv(network_path)
        prefix = args.output / "results" / name
        log_path = args.output / "reports" / f"{name}_alignment.log"
        with log_path.open("w", encoding="utf-8") as log:
            subprocess.run([args.rscript, str(MODULE / "scripts" / "compare_networks.R"),
                            str(wt_path), str(network_path), args.target, str(prefix),
                            args.r_library, str(args.ma_dim), str(args.seed)],
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        dr = pd.read_csv(f"{prefix}_dr.csv")
        delta = network.to_numpy() - wt.to_numpy()
        metrics.append({"condition": name, "changed_edges": int(np.count_nonzero(delta)),
                        "delta_frobenius": float(np.linalg.norm(delta)),
                        "median_distance": float(dr.distance.median()),
                        "significant_non_target": int(((dr["p.adj"] < .05) & (dr.gene != args.target)).sum()),
                        "tested_genes": len(dr)})
    pd.DataFrame(metrics).to_csv(args.output / "results" / "sensitivity_metrics.csv", index=False)
    metadata = {"status": "RUN_OK", "scientific_status": "Experimental",
                "target": args.target, "genes": len(wt), "folds": args.folds,
                "ma_dim": args.ma_dim, "seed": args.seed,
                "max_genes": args.max_genes,
                "ma_method": "scTenifoldNet::manifoldAlignment",
                "dr_method": "scTenifoldKnk::dRegulation; empiricalNull=FALSE; direction=NULL",
                "orientation": "row=source; column=target; no additional transposition",
                "upstream_oe_branch": "Level 1; positive target row edges only",
                "upstream_oe_revision": "860d0ea59dc0aeb564a286efae00f562ce84e138",
                "wt_sha256": hashlib.sha256(args.wt.read_bytes()).hexdigest(),
                "ko_sha256": hashlib.sha256(args.ko.read_bytes()).hexdigest(),
                "elapsed_seconds": time.perf_counter() - started,
                "versions": {p: importlib.metadata.version(p) for p in ["numpy", "pandas", "scipy", "torch", "virt_overexpression", "virt_knock"]}}
    (args.output / "reports" / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

