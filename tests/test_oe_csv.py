"""CSV identity regression without importing upstream OE or torch."""
import importlib.util
from io import StringIO
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np
import pandas as pd

PATH = Path(__file__).resolve().parents[1] / "modules/02_scTenifoldOE/scripts/run_oe.py"
SPEC = importlib.util.spec_from_file_location("run_oe_csv", PATH)
OE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(OE)


class TestNetworkCSV(unittest.TestCase):
    def test_numeric_ids_keep_leading_zeros_and_numeric_weights(self):
        text = ",001,002\n001,0,0.125\n002,0.5,0\n"
        network = OE.read_network_csv(StringIO(text))
        self.assertEqual(network.index.tolist(), ["001", "002"])
        self.assertEqual(network.columns.tolist(), ["001", "002"])
        self.assertTrue(all(pd.api.types.is_numeric_dtype(dtype) for dtype in network.dtypes))
        np.testing.assert_array_equal(network.to_numpy(), [[0, .125], [.5, 0]])
        np.testing.assert_array_equal(network.to_numpy(), pd.read_csv(StringIO(text), index_col=0).to_numpy())
        OE.validate_network(network, "001")

    def test_literal_na_id_survives_read_and_export_roundtrip(self):
        text = ",001,002,NA\n001,0,0.125,-0.25\n002,0.5,0,0.75\nNA,-1.25,2.5,0\n"
        network = OE.read_network_csv(StringIO(text))
        self.assertEqual(network.index.tolist(), ["001", "002", "NA"])
        self.assertEqual(network.columns.tolist(), ["001", "002", "NA"])
        self.assertFalse(network.index.isna().any())
        np.testing.assert_array_equal(network.to_numpy(), [[0, .125, -.25], [.5, 0, .75], [-1.25, 2.5, 0]])
        OE.validate_network(network, "NA")
        reloaded = OE.read_network_csv(StringIO(network.to_csv()))
        pd.testing.assert_frame_equal(reloaded, network)

    def test_literal_na_target_is_excluded_from_response_count(self):
        # Stub R output/identity dose only; execute real CLI reads and target exclusion.
        text = ",001,002,NA\n001,0,0.125,-0.25\n002,0.5,0,0.75\nNA,-1.25,2.5,0\n"
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            wt = OE.read_network_csv(StringIO(text))
            ko = wt.copy()
            ko.loc["NA"] = 0
            np.testing.assert_array_equal(ko.loc["NA"].to_numpy(), [0, 0, 0])
            wt.to_csv(root / "WT.csv")
            ko.to_csv(root / "KO.csv")

            def emit_dr(command, **kwargs):
                pd.DataFrame({"gene": ["001", "002", "NA"], "distance": [.1, .2, 1.],
                              "p.adj": [.01, 1., .01]}).to_csv(f"{command[5]}_dr.csv", index=False)

            argv = ["run_oe.py", "--wt", str(root / "WT.csv"), "--ko", str(root / "KO.csv"),
                    "--target", "NA", "--folds", "1", "--rscript", "unused",
                    "--r-library", "unused", "--output", str(root / "output")]
            with mock.patch("sys.argv", argv), mock.patch.object(OE.subprocess, "run", side_effect=emit_dr), \
                 mock.patch.object(OE, "apply_oe", side_effect=lambda network, *args: network.copy()), \
                 mock.patch.object(OE.importlib.metadata, "version", return_value="stub"):
                OE.main()
            metrics = pd.read_csv(root / "output/results/sensitivity_metrics.csv")
            self.assertEqual(metrics.significant_non_target.tolist(), [1, 1])


if __name__ == "__main__":
    unittest.main()
