#!/usr/bin/env python3
"""Check release language, local links, source-ID lists, and source archive contents."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tarfile


ROOT = Path(__file__).resolve().parents[1]
HAN = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff\uf900-\ufaff\U00020000-\U0002fa1f]")
REQUIRED = (
    "assets/metalprot-logo.jpg",
    "examples/basic_search/demo.ipynb",
    "data/pdb_ids/zn_pdb_ids.txt",
    "data/pdb_ids/fe_pdb_ids.txt",
    "data/pdb_ids/metadata.json",
)


def source_files():
    if (ROOT / ".git").exists():
        names = subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT).decode().split("\0")
        return [ROOT / name for name in names if name]
    excluded = {"Database", "outputs", "build", "dist", ".pytest_cache", ".venv", "__pycache__"}
    return [p for p in ROOT.rglob("*") if p.is_file()
            and not any(part in excluded or part.endswith(".egg-info")
                        for part in p.relative_to(ROOT).parts)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sdist", type=Path, help="Also check a built source distribution")
    args = parser.parse_args()
    errors = []
    files = source_files()
    for path in files:
        relative = path.relative_to(ROOT)
        name = relative.as_posix()
        if relative.parts[0] in {"Database", "legacy", "release_backups", "outputs"}:
            errors.append("Private data or archived material in public files: " + name)
        if relative.parts[0] == "data" and name not in REQUIRED:
            errors.append("Unexpected public data file: " + name)
        if HAN.search(name):
            errors.append("Non-English filename: " + name)
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        if path.suffix in {".json", ".ipynb"}:
            text = json.dumps(json.loads(text), ensure_ascii=False)
        if HAN.search(text):
            errors.append("Chinese text in release file: " + name)
        if path.suffix == ".md":
            links = re.findall(r"\]\(([^)]+)\)", text) + re.findall(r'<img[^>]+src="([^"]+)"', text)
            for link in links:
                if re.match(r"^[a-zA-Z][\w+.-]*:", link) or link.startswith("#"):
                    continue
                target = link.split("#", 1)[0]
                if not (path.parent / target).exists():
                    errors.append("Broken local link in " + name + ": " + link)
    for name in REQUIRED:
        if not (ROOT / name).is_file():
            errors.append("Missing release asset: " + name)
    metadata = json.loads((ROOT / "data/pdb_ids/metadata.json").read_text())
    for label, entry in metadata["libraries"].items():
        path = ROOT / entry["pdb_id_file"]
        data = path.read_bytes()
        ids = data.decode().splitlines()
        if len(ids) != entry["unique_pdb_id_count"] or ids != sorted(set(ids)):
            errors.append(label + " ID count, sorting, or uniqueness mismatch")
        if not all(re.fullmatch(r"[0-9][A-Z0-9]{3}", value) for value in ids):
            errors.append(label + " has an invalid PDB ID")
        if hashlib.sha256(data).hexdigest() != entry["pdb_id_file_sha256"]:
            errors.append(label + " ID list checksum mismatch")
    if args.sdist:
        with tarfile.open(args.sdist) as archive:
            names = {member.name.split("/", 1)[-1] for member in archive.getmembers() if member.isfile()}
        for name in REQUIRED:
            if name not in names:
                errors.append("Source archive is missing " + name)
        for name in names:
            if name.startswith(("Database/", "legacy/", "release_backups/", "outputs/")):
                errors.append("Unexpected source archive content: " + name)
    print(json.dumps({"checked_files": len(files), "errors": errors}, indent=2))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
