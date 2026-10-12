"""Checks of saved real-data results; no synthetic biological result is generated."""
import json
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
MODULE = REPO / "modules/02_scTenifoldOE"


@unittest.skipUnless((MODULE / "reports/run_metadata.json").exists(), "Run the real OE analysis first")
class TestOERealResults(unittest.TestCase):
    def test_null_control_and_complete_finite_results(self):
        for name in ["KO", "OE_1", "OE_1.5", "OE_2", "OE_3"]:
            table = pd.read_csv(MODULE / "results" / f"{name}_dr.csv")
            self.assertEqual(len(table), 600)
            self.assertTrue(table.gene.is_unique)
            self.assertTrue(np.isfinite(table[["distance", "FC", "p.value", "p.adj"]]).all().all())
            self.assertTrue(table["p.adj"].between(0, 1).all())
        null = pd.read_csv(MODULE / "results/OE_1_dr.csv")
        self.assertTrue((null["p.adj"] == 1).all())
        self.assertLess(null.distance.max(), 1e-12)

    def test_reused_KO_matches_official_module(self):
        official = pd.read_csv(REPO / "modules/01_scTenifoldKnk/results/single_SPI1.csv").set_index("gene")
        compare = pd.read_csv(MODULE / "results/KO_dr.csv").set_index("gene").loc[official.index]
        np.testing.assert_allclose(compare.distance, official.distance, rtol=1e-5, atol=1e-13)

    def test_enrichment_uses_network_background_and_full_family(self):
        table = pd.read_csv(MODULE / "results/pathway_enrichment.csv")
        self.assertTrue((table.background_genes == 599).all())
        self.assertTrue(table.p_adj.between(0, 1).all())
        self.assertTrue(((table.overlap == 0) & (table.p_value == 1)).any())
        for _, subset in table.groupby(["group", "library"]):
            self.assertEqual(subset.term.nunique(), len(subset))

    def test_figure_index_points_to_real_outputs(self):
        manifest = json.loads((MODULE / "figures/figure_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(len(manifest), 10)
        for record in manifest:
            self.assertIn("Experimental", record["validity"])
            for field in ["png", "pdf", "svg", "analysis_script", "plot_script"]:
                self.assertTrue((REPO / record[field]).is_file(), record[field])
            for source in record["source_tables"]:
                self.assertTrue((MODULE / source).is_file(), source)


if __name__ == "__main__":
    unittest.main()
