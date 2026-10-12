"""Retrieve the public Enrichr libraries recorded with this OE analysis."""
import argparse
import hashlib
import json
import urllib.request
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    module = Path(__file__).resolve().parents[1]
    expected = json.loads((module / "reports/pathway_provenance.json").read_text(encoding="utf-8"))
    args.output.mkdir(parents=True, exist_ok=True)
    for record in expected:
        path = args.output / (record["library"] + ".gmt")
        urllib.request.urlretrieve(record["source"], path)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != record["sha256"]:
            raise ValueError(f"{record['library']} downloaded checksum changed; re-audit annotation before comparison")
        print(record["library"], digest)


if __name__ == "__main__":
    main()
