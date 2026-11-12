#!/usr/bin/env bash
# Run the basic example in the active Python environment from any working directory.
set -euo pipefail
if (( $# > 1 )); then
    echo "Usage: bash run.sh [output_parent_directory]" >&2
    exit 2
fi
example_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd -- "$example_dir/../.." && pwd)"
output_dir="${1:-$repo_dir/outputs/basic_search}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1
python "$repo_dir/scripts/run_search.py" --config "$example_dir/config.json" --check-only
# Read the current result directory from the final JSON line in a temporary log.
run_log="$(mktemp)"
trap 'rm -f -- "$run_log"' EXIT
python "$repo_dir/scripts/run_search.py" --config "$example_dir/config.json" --outdir "$output_dir" | tee "$run_log"
result_dir="$(python -c 'import json, pathlib, sys; lines = pathlib.Path(sys.argv[1]).read_text().splitlines(); print(json.loads(lines[-1])["result_directory"])' "$run_log")"
python "$repo_dir/scripts/compare_example.py" --run-dir "$result_dir"
