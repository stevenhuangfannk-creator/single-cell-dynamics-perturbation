"""Check saved real-output scientific invariants, independently of the R implementation."""
from pathlib import Path
import csv
import json
import math
import unittest

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "modules/01_scTenifoldKnk"

def rows(name):
    with (MODULE/"results"/name).open(encoding="utf-8",newline="") as handle:
        return list(csv.DictReader(handle))

def network(name):
    with (MODULE/"objects"/name).open(encoding="utf-8",newline="") as handle:
        reader=csv.reader(handle)
        genes=next(reader)[1:]
        values={row[0]:[float(v) for v in row[1:]] for row in reader}
    return genes,values

class KnkResultTests(unittest.TestCase):
    def test_raw_input_provenance_and_target_selection(self):
        summary=json.loads((MODULE/"results/pbmc_input_summary.json").read_text(encoding="utf-8"))
        self.assertEqual(summary["selected_cells"],480)
        self.assertEqual(summary["selected_genes"],600)
        self.assertEqual(summary["duplicate_symbols_collapsed"],95)
        audit={r["gene"]:r for r in rows("candidate_expression_audit.csv")}
        self.assertGreaterEqual(float(audit["SPI1"]["detection_rate"]),0.05)
        self.assertLess(float(audit["NFKB1"]["detection_rate"]),0.05)
        self.assertNotIn("NFKB1",summary["qualified_targets"])

    def test_ko_changes_only_target_outgoing_row(self):
        genes,wt=network("WT.csv"); ko_genes,ko=network("KO.csv")
        self.assertEqual(genes,ko_genes); self.assertEqual(list(wt),genes)
        self.assertEqual(len(genes),600); self.assertEqual(len(set(genes)),600)
        self.assertGreater(sum(v!=0 for v in wt["SPI1"]),0)
        self.assertTrue(all(v==0 for v in ko["SPI1"]))
        for i,g in enumerate(genes):
            self.assertEqual(wt[g][i],0)
            self.assertTrue(all(math.isfinite(v) for v in wt[g]))
            if g!="SPI1": self.assertEqual(wt[g],ko[g])
        _,multi=network("multi_KO.csv")
        self.assertTrue(all(v==0 for target in ("SPI1","CEBPB") for v in multi[target]))
        for g in genes:
            if g not in ("SPI1","CEBPB"): self.assertEqual(wt[g],multi[g])

    def test_dr_fc_is_squared_network_distance_not_expression_fc(self):
        with (MODULE/"results/single_SPI1_MA.csv").open(encoding="utf-8",newline="") as handle:
            ma=list(csv.reader(handle))[1:]
        noise=math.sqrt(2.220446049250313e-16)*max(abs(float(v)) for r in ma for v in r[1:])
        for filename,excluded in [("single_SPI1.csv",{"SPI1"}),("multi_SPI1_CEBPB.csv",{"SPI1","CEBPB"})]:
            data=rows(filename)
            expectation=sum(float(r["distance"])**2 for r in data if r["gene"] not in excluded)/(len(data)-len(excluded))
            for row in data:
                fc=float(row["FC"])
                self.assertAlmostEqual(fc,float(row["distance"])**2/expectation,places=8)
                # Official implementation assigns p=1 below its MA floating-point noise floor.
                if float(row["p.value"])!=1:
                    self.assertAlmostEqual(float(row["p.value"]),math.erfc(math.sqrt(fc/2)),places=11)
                else:
                    if filename=="single_SPI1.csv": self.assertLessEqual(float(row["distance"]),noise*(1+1e-10))
                    # Joint-MA CSV is not part of the initial delivery; exact noise-floor
                    # verification above applies to single KO, joint p=1 remains preserved.
                self.assertTrue(0<=float(row["p.adj"])<=1)

    def test_negative_control(self):
        control=rows("negative_control_WT_vs_WT.csv")
        self.assertEqual(len(control),600)
        self.assertTrue(all(float(r["p.value"])==1 for r in control))
        self.assertTrue(all(float(r["p.adj"])==1 for r in control))

    def test_bounded_screen_has_expected_axes(self):
        path=MODULE/"results/candidate_screen_distance.csv"
        with path.open(encoding="utf-8",newline="") as handle:
            reader=csv.reader(handle); genes=next(reader)[1:]; matrix=list(reader)
        self.assertEqual(len(genes),600); self.assertEqual(len(matrix),9)
        self.assertEqual(len(set(r[0] for r in matrix)),9)
        self.assertTrue(all(math.isfinite(float(v)) and float(v)>=0 for r in matrix for v in r[1:]))

    def test_unstable_runs_preserved(self):
        stability=rows("stability_summary.csv")
        self.assertEqual(len(stability),6)
        self.assertEqual({(int(r["seed"]),int(r["cells"])) for r in stability},{(s,n) for s in (1,2,3) for n in (240,400)})
        independent=[float(r["spearman"]) for r in stability if r["run"]!="seed1_cells400"]
        self.assertLess(max(independent),0.3)
        sensitivity=rows("sensitivity1200/parameter_sensitivity_summary.csv")
        self.assertEqual(len(sensitivity),3)
        self.assertTrue(all(int(r["genes"])==1200 and int(r["networks"])==10 for r in sensitivity))
        pooled={r[""]:r for r in rows("stability_distances.csv")}
        for run in stability:
            for row in rows("stability_"+run["run"]+".csv"):
                if row["gene"]!="SPI1": self.assertAlmostEqual(float(row["distance"]),float(pooled[row["gene"]][run["run"]]),places=12)

    def test_ora_background_and_statistics(self):
        data=rows("enrichment_all_terms.csv")
        self.assertEqual(len(data),3626)
        self.assertEqual({r["universe"] for r in data},{"600"})
        self.assertTrue(any(int(r["overlap"])==0 for r in data))
        # Independently recompute one positive overlap exactly with integer combinations.
        row=next(r for r in data if int(r["overlap"])>=3)
        N,K,n,k=(int(row[c]) for c in ("universe","term_background","query_genes","overlap"))
        expected=sum(math.comb(K,x)*math.comb(N-K,n-x) for x in range(k,min(K,n)+1) if 0<=n-x<=N-K)/math.comb(N,n)
        self.assertAlmostEqual(float(row["p_value"]),expected,places=13)
        self.assertTrue(all(float(r["p_value"])<=float(r["p_adj"])+1e-12 for r in data))

    def test_figure_manifest_outputs(self):
        manifest=json.loads((MODULE/"figures/figure_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(len(manifest),11)
        for entry in manifest:
            for ext in ("png","pdf","svg"):
                path=ROOT/entry[ext]
                self.assertTrue(path.exists(),str(path)); self.assertGreater(path.stat().st_size,1000)
            self.assertIn("PARTIAL",entry["status"])

if __name__=="__main__": unittest.main()
