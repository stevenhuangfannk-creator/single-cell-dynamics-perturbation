"""对 A 异常执行独立原始表达/细胞身份对照，不运行或修改模型。"""
import argparse
import gzip
import json
from pathlib import Path
import h5py
import numpy as np

def audit(source, prepared, output):
    with h5py.File(prepared, 'r') as f:
        genes = f['var']['_index'].asstr()[:]
        cells = f['obs']['_index'].asstr()[:]
        c = f['obs']['celltype']
        labels = c['categories'].asstr()[:][c['codes'][:]]
        dense = np.zeros(tuple(f['X'].attrs['shape']), dtype=np.float32)
        vals, ix, ptr = f['X']['data'][:], f['X']['indices'][:], f['X']['indptr'][:]
        for i in range(len(cells)):
            dense[i, ix[ptr[i]:ptr[i+1]]] = vals[ptr[i]:ptr[i+1]]
    with gzip.open(source, 'rt') as f:
        original_cells = next(f).rstrip('\r\n').split('\t')[1:]
        next(f)
        original_labels = next(f).rstrip('\r\n').split('\t')[1:]
        cellmap = {c: i for i, c in enumerate(original_cells)}
        selected = np.array([cellmap[c] for c in cells])
        assert np.array_equal(labels, np.asarray(original_labels)[selected]), 'cell label mismatch'
        gmap = {g: i for i, g in enumerate(genes)}
        verified = 0
        for line in f:
            gene, values = line.split('\t', 1)
            if gene in gmap:
                original = np.fromstring(values, dtype=np.float32, sep='\t')[selected]
                assert np.array_equal(original, dense[:, gmap[gene]]), f'counts mismatch {gene}'
                verified += 1
    result = {'status': 'RUN_OK', 'cells_checked': len(cells), 'hvg_checked': verified,
              'counts_values_checked': len(cells)*verified,
              'cell_labels_exact_match': True, 'source_count_entries_exact_match': True,
              'interpretation': '排除预处理细胞标签错配及训练HVG计数错位；这不证明模型发育排序正确。'}
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    Path(output).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(result, ensure_ascii=True))

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--source', required=True)
    p.add_argument('--prepared', required=True)
    p.add_argument('--output', required=True)
    audit(**vars(p.parse_args()))
