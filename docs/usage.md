# Configuration and search outputs

This guide describes `scripts/run_search.py` and `metalprot/basic/filter.py`. Example values come from `examples/basic_search/config.json`; they are not universal recommendations for every metal or backbone. Distances are in angstroms and angles are in degrees.

## Command line and paths

- `--config`: required JSON configuration path, resolved from the current working directory.
- `--outdir`: output parent directory, default `outputs/basic_search`, resolved from the current working directory. Each search creates a separate subdirectory.
- `--check-only`: validates configuration, input-file existence, selected backbone atoms, and filter parameter names without running the search or loading database objects.

Within JSON, relative `target`, `database`, `geometry`, and `geometry_o` paths are resolved from the configuration file's directory. Adjust them when copying a configuration elsewhere. Absolute paths are also supported. JSON does not support comments; use `description` to describe a task.

To verify a separate local database manifest, run `python scripts/verify_database.py --manifest <manifest.json> --root <data-root>`. The root is the base for paths stored inside the manifest. This optional verification does not deserialize database objects.

## Configuration fields

All fields below are required unless marked optional. Missing top-level fields are not filled automatically.

| Field | Basic example | Meaning |
|---|---|---|
| `target` | `input.pdb` | Input protein PDB file |
| `database` | `../../Database/ZN_HDEC_20210113` | Directory containing the three processed files from one library |
| `positions` | `[["A",35],["A",61],["A",65]]` | At least three distinct PDB chain identifiers and integer residue numbers, not zero-based sequence indices |
| `num_contacts` | `3` | Number of coordinating residues per candidate; this entry point accepts only 3 |
| `metal_metal_dist` | `0.7` | Positive finite distance threshold for compatible candidate metal positions |
| `density_radius` | `0.6` | Positive finite radius for local support calculations; changing it affects scores |
| `allowed_aa_combinations` | `[["H","H","H"]]` | Nonempty list of three-residue compositions using H, D, E, or C; compositions are compared after sorting, not assigned position by position |
| `filters` | See below | Parameters passed to `Search_filter`; omitted filter fields use class defaults |
| `description` | Task description | Optional; does not affect search |
| `geometry` | Not set | Optional reference geometry PDB; defaults to the bundled tetrahedral resource |
| `geometry_o` | Not set | Optional open-site geometry PDB; the open-site clash switch controls whether that filter is active |

Each selected position must have unique N, CA, and C backbone atoms. The position format cannot separately specify insertion codes. Inputs with multiple models, alternate conformations, missing atoms, or duplicate chain/residue selections should be prepared before searching. Arbitrary complex multichain inputs have not been comprehensively validated.

Custom reference geometries must follow the core atom-selection and ordering conventions. File existence alone does not establish compatibility; the example does not validate an octahedral workflow. Second-shell search is disabled in this entry point and cannot be enabled with an unknown JSON field.

## Filter parameters

Defaults refer to the `Search_filter` constructor; example values include inherited defaults for fields not explicitly present in the example JSON.

| Parameter | Default | Example | Purpose |
|---|---|---|---|
| `filter_abple` | `false` | `false` | Filter by ABPLE backbone conformation class |
| `filter_phipsi` | `true` | `true` | Filter candidate backbone phi/psi angles |
| `max_phipsi_val` | `15` | `25` | Phi/psi angular tolerance |
| `filter_vdm_score` | `false` | `false` | Enable a single-motif score cutoff |
| `min_vdm_score` | `0` | `0` | Score cutoff when enabled |
| `filter_vdm_count` | `false` | `false` | Enable a cluster-size cutoff |
| `min_vdm_clu_num` | `20` | `20` | Minimum cluster size when enabled |
| `after_search_filter_geometry` | `false` | `true` | Filter geometry after candidate combination search |
| `filter_based_geometry_structure` | `false` | `false` | Use reference geometry and tolerances instead of explicit distance/angle ranges |
| `angle_tol` | `5` | `5` | Angular tolerance in reference-geometry mode |
| `aa_aa_tol` | `0.5` | `0.5` | Coordinating-atom distance tolerance in reference-geometry mode |
| `aa_metal_tol` | `0.2` | `0.2` | Metal-to-coordinating-atom distance tolerance in reference-geometry mode |
| `pair_angle_range` | `null` | `[85,140]` | Allowed coordinating-atom/metal/coordinating-atom angles |
| `pair_aa_aa_dist_range` | `null` | `[2.8,4.5]` | Allowed distances between coordinating atoms |
| `pair_metal_aa_dist_range` | `null` | `null` | Allowed metal-to-coordinating-atom distances; no range check when unset |
| `after_search_filter_qt_clash` | `false` | `false` | Check candidate/candidate and candidate/backbone clashes |
| `vdm_vdm_clash_dist` | `2.7` | `2.7` | Candidate/candidate clash threshold |
| `vdm_bb_clash_dist` | `2.2` | `2.2` | Candidate/backbone clash threshold |
| `after_search_open_site_clash` | `false` | `false` | Check open-site clashes |
| `open_site_dist` | `3.0` | `3.0` | Open-site clash threshold |
| `write_filtered_result` | `false` | `false` | Write filtered results |
| `selfcenter_filter_member_phipsi` | `true` | `true` | Apply phi/psi checks to members used for local support |

Ranges use `[minimum, maximum]`, including endpoints in range-based filtering. A threshold is effective only when its corresponding switch and mode are enabled. The entry point does not exhaustively validate numeric ranges for every filter. Reference-geometry mode also derives distance ranges from the structure and tolerances; inspect the recorded run settings.

## Using another target

1. Copy the configuration, update input and database paths, and select complete backbone positions using PDB numbering.
2. Set coordination composition, geometry, clashes, and thresholds for the research task; check the scoring reference population described in [Methods](methods.md).
3. Run `--check-only`, then search into a new output directory.
4. Inspect structures, filter status, and scores; retain the configuration and run manifest. The bundled comparator is specific to the original 5OD1 example.

## Outputs

Each run contains candidate structures, the core parameter records, a summary TSV, and `run_manifest.json`. The summary `_summary_*.tsv` contains one row per retained representative design, not every compatible candidate.

| Field | Meaning |
|---|---|
| `Wins` | PDB chain/residue positions, for example A-35_A-61_A-65 |
| `AAs` | Coordinating residue types |
| `vdMIDs` | Cluster identifiers within each residue type, not PDB IDs |
| `DensityRadius` | Radius of the local support sphere |
| `CluScore` | Combined cluster score described in Methods |
| `OverlapScore` / `OverlapScoreLn` | Spatial support score and its natural logarithm |
| `GeoRmsd` | RMSD of the implementation's geometry alignment |
| `aa_aa_dists` / `metal_aa_dists` | Coordinating-atom and metal-to-coordinating-atom distances |
| `Pair_angles` | Coordinating-atom/metal/coordinating-atom angles |
| `overlaps#` / `clu_nums` | Local supporting member counts and cluster sizes |
| `proteinPhiPsi` / `CentroidPhiPsi` | Backbone and candidate phi/psi angles |

The output schema includes second-shell fields and filter flags even when the corresponding functionality is disabled. Default values do not mean a check has passed. The example disables second-shell search, side-chain clashes, and open-site clashes.

`run_manifest.json` records resolved input paths, input hashes, configuration, library size, scoring reference values, Python/dependency versions, platform, elapsed time, and output files. It supports run provenance, not independent validation of the entire paper.
