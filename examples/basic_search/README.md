# Three-histidine search at the 5OD1 zinc site

This example checks installation, Metal Mer backbone alignment, compatibility search, scoring, and output generation. It is not the 2CBA reconstruction in Fig. 1g or a validation of experimental activity.

## Inputs

`input.pdb` is an unchanged copy of the existing 5OD1 fixture in `tests/test_data/5od1_zn.pdb`. `config.json` selects chain A residues 35, 61, and 65, allows H/H/H combinations, and uses the local `ZN_HDEC_20210113` database. `expected/` contains the reference summary and its original environment record.

Activate the `metalprot` environment created using the root README. Before running, provide `AAMetalPhiPsi.pkl`, `all_metal_vdm.pkl`, and `cluster_centroid_dict.pkl` from the original Zn library, and set `database` in `config.json`. The default `../../Database/ZN_HDEC_20210113` is a local path convention, not a download location. The example does not download data.

If the database provider supplies a SHA-256 manifest, verification is optional:

```bash
python scripts/verify_database.py --manifest Database/manifest.json --root .
```

Both the database and manifest are prepared locally. Using a different database can change motif counts, cluster identifiers, and scores, so the original reference comparison may no longer apply.

## Run

From the repository root:

```bash
bash examples/basic_search/run.sh
```

You can also invoke the script by absolute path from another directory. Supply an optional output parent directory as its only argument:

```bash
bash examples/basic_search/run.sh /tmp/metalprot-example
```

The script uses the active environment's Python and sets BLAS, OpenMP, and Numba thread counts to one. Each invocation creates its own result directory and compares only that run. The default output parent is `outputs/basic_search/` in the repository.

For interactive use, open `demo.ipynb` with the installed environment as the kernel. It finds the repository root from either the example directory or the root directory and runs the same validated entry points.

## Expected results

Input validation reports `status: inputs_valid`. With the original input, settings, and Zn library, the search produces 3,814 compatible candidates before post-search filtering and retains four representative designs. Typical search time is about 30-35 seconds on the validated CPU setup; runtime depends on the machine.

A successful comparison reports `design_count: 4`, `errors: []`, and exits with status zero. The comparator checks design keys, cluster sizes, and local support counts, with `rel_tol=1e-6` and `abs_tol=0.011` for CluScore, OverlapScore, OverlapScoreLn, and GeoRmsd. It does not compare timestamps or require byte-identical output PDBs. The reference was originally generated using NumPy 1.26.4; the supported installation pins 1.23.5 to satisfy ProDy 2.4.1.

The configuration enables phi/psi and geometric filtering, but disables ABPLE, single-motif score/count, side-chain clash, and open-site clash filters. Four retained designs do not imply that every possible design filter has been passed. Consult [Usage](../../docs/usage.md) before changing settings. Reference comparison applies only to this original example.
