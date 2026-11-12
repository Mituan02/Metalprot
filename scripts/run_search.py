#!/usr/bin/env python3
"""Run a first-shell search from JSON and record inputs, environment, and results."""

import argparse
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import platform
import sys
import time
import uuid


DATABASE_FILES = ("all_metal_vdm.pkl", "AAMetalPhiPsi.pkl", "cluster_centroid_dict.pkl")


def sha256(path):
    """Hash a file in chunks without loading the entire database into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_config(path):
    """Resolve inputs relative to the configuration; reject unknown or missing fields."""
    path = Path(path).resolve()
    config = json.loads(path.read_text(encoding="utf-8"))
    required = {"target", "database", "positions", "num_contacts", "metal_metal_dist",
                "density_radius", "allowed_aa_combinations", "filters"}
    optional = {"description", "geometry", "geometry_o"}
    if not isinstance(config, dict) or not required <= config.keys():
        raise ValueError("Configuration is missing required search fields")
    if set(config) - required - optional:
        raise ValueError("Unknown configuration fields: " + ", ".join(sorted(set(config) - required - optional)))
    for name in ("target", "database", "geometry", "geometry_o"):
        if config.get(name) is not None:
            config[name] = str((path.parent / config[name]).resolve())
    for name in ("metal_metal_dist", "density_radius"):
        if not isinstance(config[name], (int, float)) or not math.isfinite(config[name]) or config[name] <= 0:
            raise ValueError(name + " must be a finite positive number")
    if config["num_contacts"] != 3:
        raise ValueError("This entry point supports searches with three coordinating residues")
    positions = config["positions"]
    if not isinstance(positions, list) or len(positions) < 3:
        raise ValueError("Specify at least three distinct backbone positions")
    for item in positions:
        if not isinstance(item, list) or len(item) != 2 or not isinstance(item[0], str) or not isinstance(item[1], int):
            raise ValueError("Each position must be [chain identifier, PDB residue number]")
    if len({tuple(p) for p in positions}) != len(positions):
        raise ValueError("Backbone positions must be unique")
    if not isinstance(config["filters"], dict):
        raise ValueError("filters must be a parameter dictionary")
    combinations = config["allowed_aa_combinations"]
    if not isinstance(combinations, list) or not combinations:
        raise ValueError("Specify the allowed coordinating residue combinations")
    if any(len(c) != 3 or any(aa not in "HDEC" or len(aa) != 1 for aa in c) for c in combinations):
        raise ValueError("Each combination must contain three residues chosen from H, D, E, or C")
    inputs = [Path(config["target"])] + [Path(config["database"]) / n for n in DATABASE_FILES]
    inputs += [Path(config[n]) for n in ("geometry", "geometry_o") if config.get(n)]
    missing = [str(p) for p in inputs if not p.is_file()]
    if missing:
        raise FileNotFoundError("Missing inputs: " + ", ".join(missing) + ". Database files are not distributed with the code; prepare a local database and check the configured paths.")
    return config, inputs


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--outdir", type=Path, default=Path("outputs/basic_search"))
    parser.add_argument("--check-only", action="store_true", help="Check configuration, backbone positions, and database files without searching")
    args = parser.parse_args(argv)
    try:
        config, inputs = load_config(args.config)
        import prody as pr
        from metalprot.basic.filter import Search_filter
        from metalprot.search.search_selfcenter import Search_selfcenter, run_search_selfcenter

        target = pr.parsePDB(config["target"])
        if target is None or target.select("name CA") is None:
            raise ValueError("The input must contain a protein backbone")
        for chain, residue in config["positions"]:
            selected = target.select("chid " + chain + " and resnum " + str(residue) + " and name N CA C")
            if selected is None or sorted(selected.getNames()) != ["C", "CA", "N"]:
                raise ValueError(f"Position {chain}:{residue} must have unique N/CA/C backbone atoms")
        search_filter = Search_filter(**config["filters"])
        if args.check_only:
            print(json.dumps({"status": "inputs_valid", "inputs": [str(p) for p in inputs]}, ensure_ascii=False, indent=2))
            return 0

        import pickle
        started = time.perf_counter()
        database = Path(config["database"])
        objects = []
        # Load only trusted local databases supplied by the authors.
        for name in DATABASE_FILES:
            with (database / name).open("rb") as handle:
                objects.append(pickle.load(handle))
        all_metals, vdms, centroids = objects
        prefix = args.outdir.resolve() / ("search_" + uuid.uuid4().hex[:8])
        search = Search_selfcenter(
            config["target"], str(prefix), vdms, centroids, all_metals,
            num_contact_vdms=[config["num_contacts"]], metal_metal_dist=config["metal_metal_dist"],
            win_filtered=[tuple(p) for p in config["positions"]], search_filter=search_filter,
            geometry_path=config.get("geometry"), geometry_o_path=config.get("geometry_o"),
            density_radius=config["density_radius"], search_2ndshell=False,
            allowed_aa_combinations=config["allowed_aa_combinations"], output_wincomb_overlap=False,
        )
        run_search_selfcenter(search)
        result_dir = Path(search.workdir)
        report = {
            "status": "completed", "elapsed_seconds": time.perf_counter() - started,
            "retained_design_count": len(search.best_aa_comb_dict),
            "database_motif_count": len(vdms), "config": config,
            "config_sha256": sha256(args.config),
            "input_sha256": {str(p): sha256(p) for p in inputs},
            "python": sys.version, "platform": platform.platform(),
            "dependencies": {n: importlib.metadata.version(n) for n in
                ["numpy", "pandas", "scipy", "scikit-learn", "numba", "ProDy", "matplotlib"]},
            "scoring_reference": search.aa_vdm_info_dict,
            "result_files": sorted(str(p.relative_to(result_dir)) for p in result_dir.rglob("*") if p.is_file()),
            "scope": "Three-residue first-shell search; not a validation of the entire paper or experimental activity",
        }
        (result_dir / "run_manifest.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"result_directory": str(result_dir), "retained_design_count": report["retained_design_count"]}, ensure_ascii=False))
        return 0
    except (ValueError, TypeError, OSError, ImportError) as exc:
        parser.exit(2, "Search failed: " + str(exc) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
