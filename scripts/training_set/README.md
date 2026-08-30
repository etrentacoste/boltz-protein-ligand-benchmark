# Training-reference manifest construction

This directory contains the script used to construct the PLINDER-based
reference manifests employed in the Boltz-1 and Boltz-2 similarity analyses.

## Script

`build_plinder_reference_manifests.py`

The script reads a PLINDER-derived similarity table containing:

- `target_system`
- `target_release_date`

It applies the following strict temporal cutoffs:

| Reference | Temporal rule |
|---|---|
| Boltz-1 | `target_release_date < 2021-09-30` |
| Boltz-2 | `target_release_date < 2023-06-01` |

For each reference population, the script generates:

- a compressed list of PLINDER system identifiers;
- a compressed list of unique PDB identifiers;
- a CSV summary containing the cutoff, counts and release-date range.

The expected reference populations are:

| Reference | PLINDER systems | Unique PDB entries |
|---|---:|---:|
| Boltz-1 | 167,997 | 53,487 |
| Boltz-2 | 200,678 | 62,659 |

## Usage

```bash
python build_plinder_reference_manifests.py \
    --similarity-parquet /path/to/all_similarity_scores.parquet \
    --output-dir /path/to/training_set \
    --temp-dir /path/to/duckdb_temporary_directory
```

The temporary directory should be located on a filesystem with sufficient
free space because DuckDB may create intermediate files while extracting and
sorting the reference identifiers.

## Interpretation

These are PLINDER-based proxy reference universes defined using temporal
cutoffs. They are not official training manifests released by the Boltz
developers.

The generated manifests used in this study are provided in the repository's
`training_set/` directory.
