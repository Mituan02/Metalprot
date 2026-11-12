#!/usr/bin/env python3
"""Verify database files against a local manifest without unpickling them."""

import argparse
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--manifest", type=Path, required=True, help="Local SHA-256 manifest, distributed separately from the code")
    args = parser.parse_args()
    root = args.root.resolve()
    try:
        manifest = json.loads(args.manifest.read_text())
    except (OSError, ValueError) as exc:
        parser.error("Cannot read the local manifest: " + str(exc))
    errors = []
    for record in manifest["files"]:
        path = root / record["path"]
        if not path.is_file():
            errors.append(record["path"] + ": missing")
            continue
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        if path.stat().st_size != record["bytes"] or digest.hexdigest() != record["sha256"]:
            errors.append(record["path"] + ": size or checksum mismatch")
    print(json.dumps({"checked": len(manifest["files"]), "errors": errors}, ensure_ascii=False, indent=2))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
