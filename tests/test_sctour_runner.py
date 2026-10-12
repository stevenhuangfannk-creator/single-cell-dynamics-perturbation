"""已训练的真实 Smoke checkpoint 回归检查；不会重新训练模型。"""
import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
MODEL_DIR = os.environ.get('SCTOUR_SMOKE_DIR')

@unittest.skipUnless(MODEL_DIR, 'set SCTOUR_SMOKE_DIR to an executed smoke_A directory')
class TestCheckpointRestore(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import anndata as ad
        scripts = ROOT/'modules/03_scTour/scripts'
        sys.path.insert(0, str(scripts))
        spec = importlib.util.spec_from_file_location('run_sctour', scripts/'run_sctour.py')
        cls.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.module)
        cls.original = ad.read_h5ad(Path(MODEL_DIR)/'inferred.h5ad')
        cls.checkpoint = str(Path(MODEL_DIR)/'model.pth')

    def test_restore_rejects_mismatched_ids_and_counts(self):
        for change in ['ids', 'counts']:
            with self.subTest(change=change), tempfile.TemporaryDirectory() as folder:
                path = Path(folder)
                data = self.original.copy()
                if change == 'ids':
                    data.obs_names = ['other_' + name for name in data.obs_names]
                else:
                    data.X[0, 0] += 1
                data.write_h5ad(path/'input.h5ad')
                with self.assertRaisesRegex(AssertionError, 'Checkpoint'):
                    self.module.run(ROOT/'modules/03_scTour/config/paper_A.json', path/'input.h5ad',
                                    path/'output', smoke=True, load=self.checkpoint)

    def test_empty_stages_are_not_assessed(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)
            cfg = json.loads((ROOT/'modules/03_scTour/config/paper_A.json').read_text())
            cfg['stage_order'] = []
            (path/'config.json').write_text(json.dumps(cfg))
            self.original.write_h5ad(path/'input.h5ad')
            self.module.run(path/'config.json', path/'input.h5ad', path/'output',
                            smoke=True, load=self.checkpoint)
            summary = json.loads((path/'output/run_summary.json').read_text(encoding='utf-8'))
            self.assertEqual(summary['developmental_direction_assessment'], 'NOT_ASSESSED')
            self.assertIsNone(summary['developmental_median_rho_before'])
            self.assertIsNone(summary['developmental_medians_strictly_increasing'])
            self.assertFalse(summary['manual_reverse'])

if __name__ == '__main__':
    unittest.main()
