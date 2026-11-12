# Metal Mer representation and METALPROT scoring

## Scope and implementation

Metal Mer and METALPROT are described in [De novo design of metalloenzymes using Metal Mer](https://doi.org/10.64898/2026.10.13.758590).

A Metal Mer is a single-residue coordination motif extracted and clustered from known metalloprotein structures. It relates local backbone coordinates, coordinating side-chain conformation, and metal position. The representation builds on the van der Mer concept; the implementation uses names such as `VDM` and `vdms`.

METALPROT aligns candidate motifs to a target backbone, uses backbone phi/psi and metal-position compatibility to construct a graph, and searches for compatible three-residue combinations. Post-search geometry and local spatial support contribute to filtering and scoring. The paper uses Zn and Fe systems; the bundled automated example validates a three-histidine Zn search. Backbone generation, substrate modeling, Rosetta sequence design, and experimental assays belong to the wider research workflow.

## Code map

- `metalprot/basic/vdmer.py`: motif representation and coordination geometry.
- `metalprot/database/`: extraction, clustering, and database utilities.
- `metalprot/search/search.py`: backbone alignment, combination scoring, and summary output.
- `metalprot/search/find_path_by_matrix.py`: compatibility matrices and fully connected combinations.
- `metalprot/search/search_selfcenter.py`: self-centered search and metal-cloud overlap.
- `metalprot/combs/`: ligand and second-shell tools outside the validated example's scope.

## Distinct scoring quantities

The database cluster score, combination score, and overlap score use different reference quantities:

1. **Database cluster score:** `database_cluster.get_clu_info_write` uses `ln(n_cluster / mean_cluster_size)`. The grouping depends on the clustering input; different libraries need not share the same reference population.
2. **Combination CluScore:** `Search_vdM.neighbor_calc_comb_score` sums `ln(n_i / m_i)` across positions, where `m_i` is the hard-coded reference median cluster size for the residue type.
3. **OverlapScore:** for local member counts `o_i` within radius `r`, reference counts `T_i`, and median sizes `m_i`, the score is `geometric_mean(o_i * T_i / m_i) / ((4/3) * pi * r^3)`. The center is the mean metal position of the aligned candidate centroids. Member counts can also depend on phi/psi filtering.

The summary includes both OverlapScore and its natural logarithm, OverlapScoreLn. These values are not binding free energies, catalytic activities, or probabilities of success. Plots using other normalization conventions need their own recorded scripts and parameters.

The scoring reference values `T/m` are H: 2662/60, E: 829/11, D: 896/22, and C: 3211/457. They are recorded in each run manifest; they are not recalculated from the 5,853 loaded Zn records. Reference populations must be checked separately when using other libraries or Fe data.

## Interpretation of the example

The 5OD1 input is a workflow check, not the paper's 2CBA result. The database may contain the same source structure. Independent generalization assessment requires explicit exclusion of test structures and homologous data.
