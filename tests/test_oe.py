"""Input and scientific invariants for the experimental OE wrapper."""
import importlib.util
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

PATH = Path(__file__).resolve().parents[1] / "modules/02_scTenifoldOE/scripts/run_oe.py"
SPEC = importlib.util.spec_from_file_location("run_oe", PATH)
OE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(OE)


class TestOE(unittest.TestCase):
    def setUp(self):
        self.wt = pd.DataFrame([[0, .4, -.2], [.1, 0, .2], [.2, .1, 0]],
                               index=["TF", "A", "B"], columns=["TF", "A", "B"])

    def test_missing_target_is_rejected(self):
        with self.assertRaises(ValueError):
            OE.validate_network(self.wt, "ABSENT")

    def test_mismatched_axes_are_rejected(self):
        with self.assertRaises(ValueError):
            OE.validate_network(self.wt.rename(columns={"A": "C"}), "TF")

    def test_nonfinite_is_rejected(self):
        self.wt.iloc[0, 1] = np.nan
        with self.assertRaises(ValueError):
            OE.validate_network(self.wt, "TF")

    def test_fallback_is_rejected(self):
        self.wt.loc["TF"] = [0, -.1, 0]
        with self.assertRaises(ValueError):
            OE.validate_network(self.wt, "TF")

    @unittest.skipUnless(importlib.util.find_spec("virt_overexpression"), "Requires pinned optional upstream package")
    def test_actual_upstream_only_amplifies_positive_target_row(self):
        altered = OE.apply_oe(self.wt, "TF", 2)
        self.assertEqual(altered.loc["TF", "A"], .8)
        self.assertEqual(altered.loc["TF", "B"], -.2)
        np.testing.assert_array_equal(altered.loc[["A", "B"]], self.wt.loc[["A", "B"]])
        np.testing.assert_array_equal(OE.apply_oe(self.wt, "TF", 1), self.wt)
        for invalid in [np.nan, np.inf, .5]:
            with self.assertRaises(ValueError):
                OE.apply_oe(self.wt, "TF", invalid)


if __name__ == "__main__":
    unittest.main()
