"""Download fixed upstream revisions and install without redistributing sources."""
import argparse
import hashlib
import json
import subprocess
import sys
import urllib.request
import zipfile
from pathlib import Path

REVISIONS = {
    "virt_knock": "53abe39616bce36aa0fb3db8d15cc0e689cb1b54",
    "virt_overexpression": "860d0ea59dc0aeb564a286efae00f562ce84e138",
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, required=True)
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    records = []
    for name, sha in REVISIONS.items():
        url = f"https://codeload.github.com/virchow-zhang/{name}/zip/{sha}"
        archive = args.cache / f"{name}-{sha}.zip"
        urllib.request.urlretrieve(url, archive)
        digest = hashlib.sha256(archive.read_bytes()).hexdigest()
        with zipfile.ZipFile(archive) as zipped:
            zipped.extractall(args.cache)
        source = args.cache / f"{name}-{sha}"
        subprocess.run([sys.executable, "-m", "pip", "install", "--no-deps", str(source)], check=True)
        records.append({"name": name, "revision": sha, "url": url, "archive_sha256": digest})
    (args.cache / "installed_sources.json").write_text(json.dumps(records, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
