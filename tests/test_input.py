"""Synthetic fixtures test input rejection and identifier alignment only."""
import tempfile
import unittest
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite
from shared.input import read_input, validate_counts, assess_pi, export_mex


class InputContract(unittest.TestCase):
    def counts(self):
        return ad.AnnData(sparse.csr_matrix([[1, 0], [2, 3]]), obs=pd.DataFrame(index=['c2', 'c1']), var=pd.DataFrame(index=['G1', 'G2']))

    def test_reject_invalid_values(self):
        for value in [-1., np.nan, np.inf, .5]:
            with self.subTest(value=value):
                data = self.counts()
                data.X = np.array([[1., value], [2., 3.]])
                with self.assertRaises(ValueError):
                    validate_counts(data)

    def test_reject_duplicates_empty_and_missing_target(self):
        data = self.counts()
        data.var_names = ['G1', 'G1']
        with self.assertRaises(ValueError):
            validate_counts(data)
        data = self.counts()
        data.X = sparse.csr_matrix([[0, 0], [1, 1]])
        with self.assertRaises(ValueError):
            validate_counts(data)
        with self.assertRaises(ValueError):
            validate_counts(self.counts(), targets=['ABSENT'])

    def test_mex_axes_and_metadata_alignment(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            data = self.counts()
            mmwrite(path / 'matrix.mtx', data.X.T)
            (path / 'features.tsv').write_text('id1\tG1\tGene Expression\nid2\tG2\tGene Expression\n')
            (path / 'barcodes.tsv').write_text('c2\nc1\n')
            (path / 'metadata.tsv').write_text('cell\tbatch\nc1\tB1\nc2\tB2\n')
            result = read_input(path, metadata=path / 'metadata.tsv')
            np.testing.assert_array_equal(result.X.toarray(), data.X.toarray())
            self.assertEqual(list(result.obs['batch']), ['B2', 'B1'])
            self.assertEqual(validate_counts(result, targets=['G2'])['targets']['G2']['detected_cells'], 1)

    def test_h5ad_requires_explicit_counts_layer(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'input.h5ad'
            data = self.counts()
            data.write_h5ad(path)
            with self.assertRaises(ValueError):
                read_input(path)
            result = read_input(path, counts_layer='X')
            self.assertEqual(validate_counts(result)['count_sum'], 6)

    def test_audited_export_roundtrip_preserves_axes_counts_and_ids(self):
        with tempfile.TemporaryDirectory() as directory:
            source = self.counts()
            export_mex(source, directory)
            restored = read_input(directory)
            np.testing.assert_array_equal(restored.X.toarray(), source.X.toarray())
            self.assertEqual(list(restored.obs_names), ['c2', 'c1'])
            self.assertEqual(list(restored.var_names), ['G1', 'G2'])
            self.assertEqual(validate_counts(restored)['count_sum'], 6)

    def test_export_roundtrip_preserves_leading_zero_ids_and_numeric_metadata(self):
        for cells, genes in [(['002', '001'], ['001', '010']), (['002', 'NA'], ['001', 'NA'])]:
            with self.subTest(cells=cells), tempfile.TemporaryDirectory() as directory:
                path = Path(directory)
                source = self.counts()
                source.obs_names = cells
                source.var_names = genes
                source.obs['batch'] = [2, 1]
                source.obs['quality'] = [2.5, np.nan]
                export_mex(source, path)
                source.obs.iloc[::-1].to_csv(path / 'metadata.tsv', sep='\t')
                restored = read_input(path, metadata=path / 'metadata.tsv')
                self.assertEqual(list(restored.obs_names), cells)
                self.assertEqual(list(restored.var_names), genes)
                np.testing.assert_array_equal(restored.X.toarray(), source.X.toarray())
                pd.testing.assert_frame_equal(restored.obs, source.obs)
                self.assertEqual(validate_counts(restored, targets=[genes[1]])['targets'][genes[1]]['detected_cells'], 1)

    def test_10x_features_preserve_leading_zero_ids_and_symbols(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            source = self.counts()
            mmwrite(path / 'matrix.mtx', source.X.T)
            for second in ['010', 'NA']:
                (path / 'features.tsv').write_text(f'0001\t001\tGene Expression\n{second}\t{second}\tGene Expression\n')
                (path / 'barcodes.tsv').write_text(f'002\n{second}\n')
                for gene_id, expected in [('symbol', ['001', second]), ('id', ['0001', second])]:
                    with self.subTest(gene_id=gene_id, second=second):
                        restored = read_input(path, gene_id=gene_id)
                        self.assertEqual(list(restored.var_names), expected)
                        self.assertEqual(list(restored.var['gene_id']), ['0001', second])
                        self.assertEqual(list(restored.obs_names), ['002', second])
                        np.testing.assert_array_equal(restored.X.toarray(), source.X.toarray())

    def test_pi_assessment_is_not_time_inference(self):
        data = self.counts()
        for key, values in {'cell_type': ['Monocyte', 'Monocyte'], 'condition': ['Healthy', 'PI'], 'batch': ['B1', 'B2'], 'donor': ['D1', 'D2']}.items():
            data.obs[key] = values
        report = assess_pi(data, 'cell_type', 'condition', 'batch', 'donor')
        self.assertTrue(report['trajectory_evidence'].startswith('NOT_EVALUATED'))
        self.assertEqual(report['donors_by_condition'], {'Healthy': 1, 'PI': 1})


if __name__ == '__main__':
    unittest.main()
