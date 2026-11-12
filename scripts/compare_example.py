#!/usr/bin/env python3
"""Compare search results with the bundled reference, ignoring row order and timestamps."""

import argparse
import csv
import json
import math
from pathlib import Path


def read_rows(path):
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    keys = [(r["Wins"], r["AAs"], r["vdMIDs"]) for r in rows]
    if len(keys) != len(set(keys)):
        raise ValueError("Results contain duplicate design keys")
    return dict(zip(keys, rows))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    reference = read_rows(root / "examples/basic_search/expected/summary.tsv")
    files = list(args.run_dir.glob("_summary_*.tsv"))
    if len(files) != 1:
        parser.error("The run directory must contain exactly one summary TSV")
    observed = read_rows(files[0])
    errors = []
    if observed.keys() != reference.keys():
        errors.append("Retained design combinations differ from the reference")
    for key in reference.keys() & observed.keys():
        for field in ["CluScore", "OverlapScore", "OverlapScoreLn", "GeoRmsd"]:
            if not math.isclose(float(reference[key][field]), float(observed[key][field]), rel_tol=1e-6, abs_tol=0.011):
                errors.append(str(key) + ": " + field + " exceeds the comparison tolerance")
        for field in ["overlaps#", "clu_nums"]:
            if reference[key][field] != observed[key][field]:
                errors.append(str(key) + ": " + field + " differs from the reference")
    print(json.dumps({"design_count": len(observed), "errors": errors}, ensure_ascii=False, indent=2))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
