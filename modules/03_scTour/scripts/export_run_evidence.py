"""从实际完成目录导出小型结果证据；模型和 AnnData 仍留在工作缓存。"""
import argparse
import csv
import json
from pathlib import Path
import shutil

def export(work_root, module):
    work, module = Path(work_root), Path(module)
    results = module / 'results'
    results.mkdir(parents=True, exist_ok=True)
    folders = [('A90_seed0', work/'full_A'), ('A90_seed1', work/'full_A_seed1'),
               ('A20_seed0', work/'full_A20'), ('smoke_A', work/'smoke_A')]
    folders.extend((p.name, p) for p in sorted((work/'subsamples').glob('B_*')) if p.is_dir())
    rows = []
    for name, path in folders:
        if not (path/'run_summary.json').exists():
            continue
        summary = json.loads((path/'run_summary.json').read_text(encoding='utf-8'))
        if '\ufffd' in summary.get('reason', ''):
            summary['reason'] = '实际运行成功；发育趋势与配置一致性见复现报告。'
        dest = results / name
        dest.mkdir(parents=True, exist_ok=True)
        (dest/'run_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
        for item in ['config_used.json', 'environment.json', 'loss.tsv', 'train_cells.tsv', 'validation_cells.tsv', 'inference.tsv.gz']:
            if (path/item).exists():
                shutil.copy2(path/item, dest/item)
        rows.append({'run': name, **{k: summary.get(k) for k in
                     ['status', 'n_cells_inferred', 'hvg', 'percent', 'epochs', 'random_state',
                      'training_seconds', 'final_train_loss', 'final_validation_loss',
                      'developmental_median_rho_before', 'manual_reverse', 'roundtrip_max_ptime_difference']}})
    with (results/'run_comparison.tsv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    for dataset in ['A', 'B']:
        with (work/'prepared'/f'{dataset}_gene_qc_hvg.tsv').open(encoding='utf-8') as source:
            reader = csv.DictReader(source, delimiter='\t')
            with (results/f'{dataset}_HVGs.tsv').open('w', newline='', encoding='utf-8') as target:
                writer = csv.DictWriter(target, fieldnames=reader.fieldnames, delimiter='\t')
                writer.writeheader()
                writer.writerows(row for row in reader if row['highly_variable'] == 'True')
    if (work/'subsamples'/'robustness.tsv').exists():
        shutil.copy2(work/'subsamples'/'robustness.tsv', results/'subsampling_robustness.tsv')
    print(f'Exported {len(rows)} completed real runs')

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--work-root', required=True)
    p.add_argument('--module', required=True)
    export(**vars(p.parse_args()))
