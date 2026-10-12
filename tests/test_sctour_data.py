"""科研输入约束与实际输出验收；在 scTour 环境运行。"""
import importlib.util
import json
import os
from pathlib import Path
import unittest

try:
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy
    HAVE_RUNTIME = True
except ImportError:
    HAVE_RUNTIME = False

ROOT = Path(__file__).resolve().parents[1]

@unittest.skipUnless(HAVE_RUNTIME, 'requires the documented scTour environment')
class TestCounts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location('prepare_geo', ROOT/'modules/03_scTour/scripts/prepare_geo.py')
        cls.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.module)

    def matrix(self):
        return ad.AnnData(np.ones((30, 30), dtype=np.float32),
                          obs=pd.DataFrame(index=[f'cell{i}' for i in range(30)]),
                          var=pd.DataFrame(index=[f'gene{i}' for i in range(30)]))

    def test_integer_counts_accepted(self):
        self.module.validate_counts(self.matrix())

    def test_invalid_input_rejected(self):
        for value in [-1, float('nan'), float('inf'), .5]:
            a = self.matrix()
            a.X[0, 0] = value
            with self.assertRaises(AssertionError):
                self.module.validate_counts(a)
        a = self.matrix()
        a.obs_names = ['duplicate'] * a.n_obs
        with self.assertRaises(AssertionError):
            self.module.validate_counts(a)

    @unittest.skipUnless(os.environ.get('SCTOUR_RESULT_DIR'), 'set SCTOUR_RESULT_DIR for executed-object validation')
    def test_actual_saved_inference(self):
        folder = Path(os.environ['SCTOUR_RESULT_DIR'])
        summary = json.loads((folder/'run_summary.json').read_text(encoding='utf-8'))
        a = ad.read_h5ad(folder/'inferred.h5ad')
        self.module.validate_counts(a)
        self.assertEqual(a.n_obs, summary['n_cells_inferred'])
        self.assertTrue(a.obs.ptime.between(0, 1).all())
        self.assertEqual(summary['roundtrip_max_ptime_difference'], 0)
        for rep in ['X_TNODE', 'X_VF']:
            self.assertEqual(a.obsm[rep].shape, (a.n_obs, 5))
            self.assertTrue(np.isfinite(a.obsm[rep]).all())
        self.assertEqual(a.obsm['X_umap_pca'].shape, (a.n_obs, 2))
        self.assertEqual(a.obsm['X_umap_latent'].shape, (a.n_obs, 2))

if __name__ == '__main__':
    unittest.main()
