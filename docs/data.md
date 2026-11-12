# Metal Mer data

Metal Mer libraries contain coordination motifs for metals including Zn, Fe, Co, Mn, and Cu. A library provides the structural motifs and cluster information used during a search.

## Source PDB IDs

The following lists contain the unique source PDB identifiers represented in the final curated libraries:

| Library | Metal Mer records | Unique PDB IDs | List |
|---|---:|---:|---|
| `ZN_HDEC_20210113` | 5,853 | 1,366 | [Zn IDs](../data/pdb_ids/zn_pdb_ids.txt) |
| `FE_HDE_20220116` | 976 | 223 | [Fe IDs](../data/pdb_ids/fe_pdb_ids.txt) |

The Fe list is derived from the Fe-only library. IDs are parsed from the query titles in `AAMetalPhiPsi.pkl`, uppercased, deduplicated, and sorted. They describe the final libraries rather than all structures initially downloaded for database construction. Library and list checksums are provided in [metadata.json](../data/pdb_ids/metadata.json).

## Preparing a search database

Database binaries are not included in this repository. Obtain the processed Metal Mer library separately and set the `database` field in the search configuration to a local directory containing all three files from the same library:

```text
AAMetalPhiPsi.pkl
all_metal_vdm.pkl
cluster_centroid_dict.pkl
```

The default example path is `Database/ZN_HDEC_20210113` relative to the repository root. You can use another local path by updating the configuration. The example does not download a database; the source PDB ID lists alone are not sufficient to run a search.

If your data provider supplies a SHA-256 manifest, verify the files with:

```bash
python scripts/verify_database.py --manifest /path/to/manifest.json --root /path/to/data-root
```

`--root` specifies the base directory for relative paths stored in the manifest. The search also records input hashes and library details in its run manifest. Load only databases from a trusted source.

The code's MIT license does not automatically apply to separately supplied databases or third-party data.
