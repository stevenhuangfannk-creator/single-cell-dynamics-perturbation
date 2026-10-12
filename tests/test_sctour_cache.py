"""旧无hash缓存的元数据回归检查，使用小型HDF5测试输入，不训练模型。"""
import importlib.util
from pathlib import Path
import tempfile
import unittest

try:
    import h5py
    import numpy as np
    HAVE_RUNTIME = True
except ImportError:
    HAVE_RUNTIME = False

ROOT = Path(__file__).resolve().parents[1]

@unittest.skipUnless(HAVE_RUNTIME, 'requires h5py and numpy')
class TestOldCache(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location('subsampling', ROOT/'modules/03_scTour/scripts/run_subsampling.py')
        cls.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.module)

    def write(self, path):
        with h5py.File(path, 'w') as f:
            for axis in ['obs', 'var']:
                group = f.create_group(axis)
                group.attrs['_index'] = '_index'
                group.create_dataset('_index', data=np.array([b'id1', b'id2']))
            for key in ['celltype', 'batch']:
                group = f['obs'].create_group(key)
                group.create_dataset('categories', data=np.array([b'one', b'two']))
                group.create_dataset('codes', data=[0, 1])
            f['obs'].create_dataset('n_genes_by_counts', data=[1, 1])
            group = f.create_group('X')
            group.create_dataset('data', data=[1., 2.])
            group.create_dataset('indices', data=[0, 1])
            group.create_dataset('indptr', data=[0, 1, 2])

    def test_same_counts_changed_annotations_are_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            source, cached = Path(folder)/'source.h5ad', Path(folder)/'cached.h5ad'
            self.write(source)
            for key in ['celltype', 'batch', 'n_genes_by_counts']:
                with self.subTest(key=key):
                    self.write(cached)
                    with h5py.File(cached, 'r+') as f:
                        if key == 'n_genes_by_counts':
                            f['obs'][key][0] = 2
                        else:
                            f['obs'][key]['codes'][0] = 1
                    with self.assertRaisesRegex(AssertionError, 'annotation/QC'):
                        self.module.verify_old_cached_input(source, cached)
            self.write(cached)
            self.module.verify_old_cached_input(source, cached)

if __name__ == '__main__':
    unittest.main()
