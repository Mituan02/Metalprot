# Testing METALPROT

## Regression tests

After installing the package, run:

```bash
python -m pytest tests/release -q
```

The 17 tests cover module imports, packaged resources, graph compatibility, native contact atoms, scoring, reference geometry, configuration validation, path resolution, and missing-database diagnostics. They do not require a Metal Mer database. `test_data/5od1_zn.pdb` supplies the small structural fixture.

GitHub Actions checks dependency consistency, package installation, regression tests, public-file requirements, and source-distribution contents. To check the public files locally, run `python scripts/check_release.py`.

## Example integration check

Prepare the original local Zn library as described in [Data](../docs/data.md), then run from the repository root:

```bash
bash examples/basic_search/run.sh
```

With the supplied input and settings, the example produces 3,814 compatible candidates before post-search filtering and retains four representative designs. A successful comparison reports `design_count: 4`, `errors: []`, and exits with status zero.

The comparison checks design keys, cluster sizes, and supporting-member counts. It uses `rel_tol=1e-6` and `abs_tol=0.011` for CluScore, OverlapScore, OverlapScoreLn, and GeoRmsd, and ignores timestamps and row order. Reference values and their generation environment are recorded in `examples/basic_search/expected/`.

The runner fixes BLAS, OpenMP, and Numba thread counts to one and compares only the current run. Public CI does not download databases or execute this full integration search.

This example validates the 5OD1 workflow, not the paper's 2CBA reconstruction or the complete Fe and second-shell design pipelines. Because the source structure may be represented in the database, it is not an independent generalization benchmark. See [Methods](../docs/methods.md) for scoring and interpretation.
