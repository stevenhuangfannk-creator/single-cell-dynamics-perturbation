"""Local ORA against audited Enrichr GMT; measured-network background is mandatory."""
from pathlib import Path
import argparse
import hashlib
import json
import urllib.request
import pandas as pd
from scipy.stats import hypergeom, false_discovery_control

LIBRARIES = ["GO_Biological_Process_2025", "KEGG_2021_Human", "Reactome_2022"]

def read_gmt(path):
    terms = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) >= 3:
            terms.setdefault(fields[0], set()).update(g.split(",")[0].strip() for g in fields[2:] if g.strip())
    return terms

def ora(terms, selected, background):
    selected = set(selected) & set(background)
    records = []
    for term, genes in terms.items():
        mapped = genes & background
        if not 3 <= len(mapped) <= 300:
            continue
        overlap = mapped & selected
        records.append(dict(term=term, overlap=len(overlap), term_background=len(mapped), query_genes=len(selected), universe=len(background), p_value=float(hypergeom.sf(len(overlap)-1, len(background), len(mapped), len(selected))), genes=";".join(sorted(overlap))))
    df = pd.DataFrame(records)
    if not df.empty:
        df["p_adj"] = false_discovery_control(df.p_value.to_numpy(), method="bh")
        df = df.sort_values(["p_adj", "p_value", "term"])
    return df

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--module", type=Path, default=Path("modules/01_scTenifoldKnk"))
    ap.add_argument("--gmt", type=Path, default=Path("work/knk/data/gmt"))
    args = ap.parse_args()
    args.gmt.mkdir(parents=True, exist_ok=True)
    background = set(pd.read_csv(args.module/"objects/pbmc_cd14_genes.tsv", header=None)[0])
    manifest = []
    output = []
    for library in LIBRARIES:
        url = f"https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={library}"
        path = args.gmt / f"{library}.gmt"
        if not path.exists():
            urllib.request.urlretrieve(url, path)
        terms = read_gmt(path)
        if len(terms) < 100:
            raise ValueError(f"Invalid GMT download: {library}")
        manifest.append(dict(library=library, url=url, bytes=path.stat().st_size, sha256=hashlib.sha256(path.read_bytes()).hexdigest(), terms=len(terms)))
        for contrast, excluded in [("single_SPI1", {"SPI1"}), ("multi_SPI1_CEBPB", {"SPI1", "CEBPB"})]:
            dr = pd.read_csv(args.module/"results"/f"{contrast}.csv")
            selected = set(dr.loc[dr["p.adj"]<0.05,"gene"]) - excluded
            df = ora(terms, selected, background)
            df["library"] = library
            df["contrast"] = contrast
            output.append(df)
    all_results = pd.concat(output, ignore_index=True)
    all_results.to_csv(args.module/"results/enrichment_all_terms.csv", index=False)
    all_results[all_results.p_adj<0.05].to_csv(args.module/"results/enrichment_significant.csv", index=False)
    summary = dict(status="RUN_OK", method="One-sided hypergeometric ORA; BH separately per library and contrast; all eligible terms including zero overlap; query=BH<0.05 non-target genes; measured background=600 network genes", background_genes=len(background), term_size_min=3, term_size_max=300, libraries=manifest, evaluated_rows=len(all_results), significant_rows=int((all_results.p_adj<0.05).sum()))
    (args.module/"results/enrichment_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))

if __name__ == "__main__":
    main()
