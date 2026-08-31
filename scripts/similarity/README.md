# Similarity calculation scripts

This directory contains the scripts used to compare the benchmark complexes
against the PLINDER-based Boltz-1 and Boltz-2 reference universes.

The processed results produced by this workflow are available under
[`../../similarity/`](../../similarity/). The reference-system manifests are
available under [`../../training_set/`](../../training_set/).

## Reference universes

| Reference | Temporal rule | PLINDER systems | Unique PDB entries |
|---|---:|---:|---:|
| Boltz-1 | release date `< 2021-09-30` | 167,997 | 53,487 |
| Boltz-2 | release date `< 2023-06-01` | 200,678 | 62,659 |

These are PLINDER-based proxy reference universes and not official training
manifests released by the Boltz developers.

## Similarity dimensions

The workflow generated six complementary similarity measurements:

1. complete-protein sequence similarity;
2. complete-protein structural similarity;
3. binding-pocket similarity;
4. protein-ligand interaction similarity;
5. three-dimensional ligand-pose similarity;
6. two-dimensional ligand chemical similarity.

Foldseek and MMseqs2 results were used as structural and sequence inputs.
PLINDER-derived scoring was used for protein, pocket and interaction
comparisons. SuCOS was used for three-dimensional ligand comparison and
Morgan fingerprints with the Tanimoto coefficient were used for chemical
similarity.

## Python scripts

### `build_boltz2_training_manifests.py`

Constructs the Boltz-2 PLINDER reference manifests using the strict temporal
cutoff:

```text
target release date < 2023-06-01
```

### `summarize_multiquery_target_manifests.py`

Combines the eligible target-system manifests generated for the benchmark
queries and reports the resulting reference-system coverage.

### `run_plinder_foldseek_scoring.py`

Combines query-system manifests, target-system manifests, Foldseek
alignments and MMseqs2 alignments to calculate PLINDER-derived protein,
pocket and interaction scores.

### `plinder_similarity_reporting.py`

Contains shared functions used when summarising PLINDER-derived similarity
scores and constructing report tables.

### `calculate_morgan_array.py`

Calculates ligand chemical similarity using Morgan fingerprints and the
Tanimoto coefficient.

### `calculate_sucos_array.py`

Calculates three-dimensional ligand similarity using SuCOS shape and
feature-overlap scores.

### `calculate_validate_sucos_shape.py`

Implements the SuCOS molecular alignment, shape-overlap and feature-overlap
functions used by `calculate_sucos_array.py`. Both scripts must remain in the
same directory because the array script imports this module directly.

### `consolidate_similarity_report.py`

Combines PLINDER, Morgan and SuCOS results into the per-complex files:

- `pairwise_scores.parquet`;
- `pdb_target_scores.parquet`;
- `similarity_summary.csv`;
- `similarity_summary.parquet`;
- `similarity_report.txt`.

## Slurm scripts

The `*_iridis.sbatch` files reproduce the resource requests and array-job
organisation used on the Iridis HPC cluster:

- `run_plinder_scoring_array_iridis.sbatch`;
- `calculate_morgan_array_iridis.sbatch`;
- `calculate_sucos_array_iridis.sbatch`;
- `consolidate_similarity_reports_iridis.sbatch`;
- `summarize_boltz2_targets_iridis.sbatch`;
- `consolidate_boltz2_reports_iridis.sbatch`.

The published scripts contain configurable paths rather than the original
user-specific Iridis paths.

## Configuration

Before submitting the Slurm jobs, define the project directory and,
optionally, the Apptainer executable.

The following is an example only. Replace both paths with real paths on the
target system:

```bash
export PLINDER_SIMILARITY_PROJECT="/absolute/path/to/plinder_similarity_project"
export APPTAINER_BIN="/absolute/path/to/apptainer"

mkdir -p "${PLINDER_SIMILARITY_PROJECT}/logs"
cd "${PLINDER_SIMILARITY_PROJECT}"
```

If `apptainer` is already available through `PATH`, the second variable may
be omitted:

```bash
export PLINDER_SIMILARITY_PROJECT="/absolute/path/to/plinder_similarity_project"
mkdir -p "${PLINDER_SIMILARITY_PROJECT}/logs"
cd "${PLINDER_SIMILARITY_PROJECT}"
```

The project directory is expected to contain the prepared PLINDER release,
query systems, reference manifests, Foldseek and MMseqs2 result tables, and
the directory hierarchy referenced by the Slurm scripts.

## Array-job submission

The array size is determined from the number of query PDB identifiers:

```bash
manifest="${PLINDER_SIMILARITY_PROJECT}/results/system_enumeration/calculable_queries/calculable_query_pdb_ids.txt"
n_queries="$(wc -l < "${manifest}")"
```

Example submission:

```bash
sbatch \
    --array="1-${n_queries}" \
    scripts/similarity/run_plinder_scoring_array_iridis.sbatch
```

Morgan and SuCOS calculations can subsequently be submitted using the same
array range:

```bash
sbatch \
    --array="1-${n_queries}" \
    scripts/similarity/calculate_morgan_array_iridis.sbatch

sbatch \
    --array="1-${n_queries}" \
    scripts/similarity/calculate_sucos_array_iridis.sbatch
```

After the component scores have completed, the final Boltz-1 reports can be
consolidated with:

```bash
sbatch \
    --array="1-${n_queries}" \
    scripts/similarity/consolidate_similarity_reports_iridis.sbatch
```

The corresponding Boltz-2 reference preparation and report consolidation
are represented by:

```bash
sbatch scripts/similarity/summarize_boltz2_targets_iridis.sbatch

sbatch \
    --array="1-${n_queries}" \
    scripts/similarity/consolidate_boltz2_reports_iridis.sbatch
```

On an actual cluster, dependencies should be added so that consolidation
starts only after all required PLINDER, Morgan and SuCOS jobs have completed.

## Data not included

The repository does not distribute:

- the complete PLINDER dataset;
- reconstructed PLINDER structures;
- Foldseek or MMseqs2 databases;
- container root filesystems;
- intermediate working directories.

These components are substantially larger than the processed outputs and
can be reconstructed from the documented PLINDER release and scripts.
