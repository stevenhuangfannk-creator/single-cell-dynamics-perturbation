"""绘图完成边界的软件回归；临时小表格不是训练或生物学证据。"""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]

try:
    import pandas as pd
    import scanpy
    import sctour
    HAVE_RUNTIME = True
except ImportError:
    HAVE_RUNTIME = False

@unittest.skipUnless(HAVE_RUNTIME, 'requires the documented scTour environment')
class TestFigureBoundaries(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location('plot_figure2', ROOT/'modules/03_scTour/scripts/plot_figure2.py')
        cls.plotter = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.plotter)

    def fake_run(self, root, percent, cells=('001', '002'), declared_percent=None):
        folder = root/f'B_{percent}'
        folder.mkdir()
        pd.DataFrame({'cell_id': cells, 'celltype': ['nIPC', 'GC-adult'],
                      'ptime': [.1, .9]}).to_csv(folder/'inference.tsv.gz', sep='\t', index=False)
        summary = {'epochs': 400, 'percent': percent/100 if declared_percent is None else declared_percent,
                   'n_cells_inferred': 2}
        (folder/'run_summary.json').write_text(json.dumps(summary), encoding='utf-8')
        return folder

    def invoke(self, root, panel):
        self.plotter.plot(None, None, root/'figures', ROOT/'modules/03_scTour/config/plot.json',
                          panel=panel, subsample_dir=root)

    def test_partial_grid_cannot_produce_h(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.fake_run(root, 95)
            with self.assertRaises(ValueError):
                self.invoke(root, 'h')
            self.assertFalse((root/'figures/figure_manifest.json').exists())

    def test_subsampling_cell_order_and_fraction_checked(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.fake_run(root, 95)
            self.fake_run(root, 20, cells=('002', '001'))
            with self.assertRaisesRegex(AssertionError, 'IDs/order'):
                self.invoke(root, 'h')
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.fake_run(root, 20, declared_percent=.95)
            with self.assertRaises(AssertionError):
                self.invoke(root, 'h')

    def test_incomplete_i_cannot_be_plotted(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root/'B_20').mkdir()
            (root/'B_20/inferred.h5ad').touch()
            with self.assertRaises(ValueError):
                self.invoke(root, 'i')

    def test_h_has_exactly_prescribed_eight_rows(self):
        frame = pd.DataFrame({'celltype': ['nIPC', 'GC-adult'], 'ptime': [.1, .9]})
        runs = {percent: frame for percent in self.plotter.SUBSAMPLING_PERCENTAGES | {30}}
        runs[20] = frame.assign(ptime=[.9, .1])
        fig, ax = self.plotter.plt.subplots()
        try:
            self.plotter.robustness(ax, runs, {'font_size': 7})
            self.assertEqual(ax.images[0].get_array().shape, (8, 2, 3))
            self.assertEqual([text.get_text() for text in ax.get_yticklabels()],
                             ['1%', '5%', '10%', '20%', '40%', '60%', '80%', '95%'])
            self.assertEqual([text.get_text() for text in ax.texts[1:]],
                             ['1.000', '1.000', '1.000', '-1.000', '1.000', '1.000', '1.000', '1.000'])
            self.assertEqual([text.get_position()[1] for text in ax.texts[1:]], list(range(8)))
        finally:
            self.plotter.plt.close(fig)

if __name__ == '__main__':
    unittest.main()
