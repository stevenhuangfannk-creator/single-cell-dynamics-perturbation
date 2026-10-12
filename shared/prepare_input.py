"""Validate and export an explicit counts object; no model fitting."""
import argparse
from pathlib import Path
import sys
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared.input import read_input, validate_counts, assess_pi, save_json

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--config', required=True)
parser.add_argument('--output', required=True)
parser.add_argument('--assess-pi', action='store_true')
args = parser.parse_args()
config = yaml.safe_load(Path(args.config).read_text(encoding='utf-8'))
adata = read_input(config['input'], config.get('counts_layer', 'counts'), config.get('metadata'), config.get('gene_id', 'symbol'))
report = validate_counts(adata, min_cells=config.get('min_cells', 1), targets=config.get('targets', []))
if args.assess_pi:
    report = assess_pi(adata, config['cell_type_key'], config['condition_key'], config['batch_key'], config['donor_key'], config.get('targets', []))
output = Path(args.output)
output.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(output / 'validated_counts.h5ad')
save_json(report, output / 'input_audit.json')
print(f'Validated {adata.n_obs} cells × {adata.n_vars} genes; audit saved to {output}')
