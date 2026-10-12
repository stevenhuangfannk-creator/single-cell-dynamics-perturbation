"""Background-aware KO/OE ORA; retain every tested pathway including zero hits."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom, spearmanr


def bh(pvals):
    pvals = np.asarray(pvals, float)
    order = pvals.argsort()
    adjusted = np.minimum.accumulate((pvals[order] * len(pvals) / np.arange(1, len(pvals) + 1))[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--module", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--gmt", type=Path, required=True)
    args = parser.parse_args()
    target = json.loads((args.module / "reports/run_metadata.json").read_text())["target"]
    tables = {key: pd.read_csv(args.module / "results" / f"{key}_dr.csv").set_index("gene") for key in ["KO", "OE_1.5", "OE_2", "OE_3"]}
    universe = set(tables["KO"].index) - {target}
    for table in tables.values():
        assert set(table.index) - {target} == universe
    sig = {k: set(t.index[t["p.adj"] < .05]) - {target} for k, t in tables.items()}
    groups = {"KO": sig["KO"], "OE 2": sig["OE_2"], "Shared": sig["KO"] & sig["OE_2"],
              "KO specific": sig["KO"] - sig["OE_2"], "OE specific": sig["OE_2"] - sig["KO"]}
    records = []
    provenance = []
    for path in sorted(args.gmt.glob("*.gmt")):
        library = path.stem
        sets = {}
        for line in path.read_text(encoding="utf-8").splitlines():
            parts = line.split("\t")
            members = {g.split(",")[0] for g in parts[2:] if g} & universe
            if 3 <= len(members) <= 300:
                sets[parts[0]] = members
        provenance.append({"library": library, "source": "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=" + library,
                           "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "tested_pathways": len(sets)})
        for group, hit_genes in groups.items():
            these = []
            for term, members in sets.items():
                overlap = members & hit_genes
                p = hypergeom.sf(len(overlap) - 1, len(universe), len(members), len(hit_genes))
                these.append({"library": library, "group": group, "term": term,
                              "background_genes": len(universe), "response_genes": len(hit_genes),
                              "pathway_genes_in_background": len(members), "overlap": len(overlap),
                              "overlap_genes": ";".join(sorted(overlap)), "p_value": float(p)})
            if these:
                for record, q in zip(these, bh([r["p_value"] for r in these])):
                    record["p_adj"] = float(q)
                    records.append(record)
    if not records:
        raise ValueError("No compatible GMT pathways; cannot report enrichment")
    result = pd.DataFrame(records)
    result.to_csv(args.module / "results/pathway_enrichment.csv", index=False)
    (args.module / "reports/pathway_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    stability = []
    reference = tables["OE_2"].loc[sorted(universe)]
    for key in ["OE_1.5", "OE_3"]:
        test = tables[key].loc[reference.index]
        top_ref, top_test = set(reference.distance.nlargest(20).index), set(test.distance.nlargest(20).index)
        stability.append({"comparison": key + " versus OE_2", "spearman_distance": spearmanr(reference.distance, test.distance).statistic,
                          "top20_jaccard": len(top_ref & top_test) / len(top_ref | top_test),
                          "scope": "Dose sensitivity on one shared WT; not across-seed biological validation"})
    pd.DataFrame(stability).to_csv(args.module / "results/dose_stability.csv", index=False)


if __name__ == "__main__":
    main()
