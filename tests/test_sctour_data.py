"""科研输入约束与实际输出验收；在 scTour 环境运行。"""
import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile
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

    def test_mex_string_ids_and_metadata_alignment(self):
        from scipy.io import mmwrite
        from scipy import sparse
        sys.path.insert(0, str(ROOT/'modules/03_scTour/scripts'))
        try:
            import prepare_custom
            rng = np.random.default_rng(17)
            counts = rng.poisson(rng.uniform(2, 9, 200), size=(60, 200))
            cells = [f'{i:03d}' for i in range(1, 61)]
            genes = ['NA'] + [f'{i:03d}' for i in range(1, 200)]
            with tempfile.TemporaryDirectory() as temporary:
                folder = Path(temporary)
                mmwrite(folder/'matrix.mtx', sparse.csr_matrix(counts.T))
                (folder/'barcodes.tsv').write_text('\n'.join(cells)+'\n', encoding='utf-8')
                (folder/'genes.tsv').write_text('\n'.join(genes)+'\n', encoding='utf-8')
                pd.DataFrame({'cell_id': cells[::-1], 'cell_type': 'Fibroblast',
                              'batch': cells[::-1]}).to_csv(folder/'metadata.tsv', sep='\t', index=False)
                prepare_custom.prepare(folder, folder/'prepared', 'cell_type', 'batch', ['Fibroblast'], 50)
                result = ad.read_h5ad(folder/'prepared/qc_counts.h5ad')
                self.assertEqual(result.obs_names.tolist(), cells)
                self.assertEqual(result.var_names.tolist(), genes)
                self.assertEqual(result.obs.batch.tolist(), cells)
                np.testing.assert_array_equal(result.X.toarray(), counts)
                self.assertEqual(result.X.dtype, np.float32)
                self.assertEqual(ad.read_h5ad(folder/'prepared/train_counts.h5ad').n_vars, 50)
        finally:
            sys.path.pop(0)

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
