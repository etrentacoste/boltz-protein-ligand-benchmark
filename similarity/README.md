# Similarity analysis

This directory contains the similarity results obtained by comparing the
100-complex benchmark against two PLINDER-based reference universes
representing structures released before the temporal cutoffs used for
Boltz-1 and Boltz-2.

These reference universes are reproducible PLINDER-based proxies. They
should not be interpreted as official training manifests released by the
Boltz developers.

## Reference universes

| Reference | Temporal rule | PLINDER systems | Unique PDB entries |
|---|---:|---:|---:|
| Boltz-1 | release date `< 2021-09-30` | 167,997 | 53,487 |
| Boltz-2 | release date `< 2023-06-01` | 200,678 | 62,659 |

The complete reference manifests and their construction procedure are
documented in [`../training_set/`](../training_set/).

## Similarity metrics

Six complementary similarity dimensions were considered:

1. **Protein sequence similarity**: similarity between the complete query
   and reference protein sequences.
2. **Protein structure similarity**: structural similarity between the
   complete query and reference proteins.
3. **Pocket similarity**: similarity between the residues forming the
   ligand-binding pockets.
4. **Protein–ligand interaction similarity**: similarity between the
   interaction environments of the query and reference ligands.
5. **Ligand pose similarity**: three-dimensional ligand similarity based
   primarily on SuCOS shape comparison.
6. **Ligand chemistry similarity**: two-dimensional chemical similarity
   calculated using Morgan fingerprints and the Tanimoto coefficient.

Foldseek was used for structural candidate generation and residue
correspondence. PLINDER-derived scores were then used for the protein,
pocket and interaction comparisons. Ligand shape and chemical similarity
were calculated separately.

## Coverage

The benchmark contains 100 protein–ligand complexes.

Complete per-target similarity reports were produced for:

| Reference | Benchmark complexes with reports | Benchmark complexes without reports |
|---|---:|---:|
| Boltz-1 | 90 | 10 |
| Boltz-2 | 90 | 10 |

The missing ten complexes did not produce a final eligible PLINDER query
system and therefore do not have per-target similarity reports. They remain
part of the structural prediction benchmark and are explicitly recorded in
`benchmark_similarity_status.csv`.

### Metric-level coverage

The 90 per-target reports for each reference set do not all contain values for
all six similarity metrics. `22MJ`, `9HZ3` and `9VQ5` have no Pocket or
Interaction similarity values. `9PSQ` has no Interaction similarity value.

These missing values reflect unavailable metric-specific comparisons, not
missing per-target report directories. Downstream analyses retain only the
available values for the relevant metric.

Absence of a report must not be interpreted as zero similarity.

## Directory structure

```text
similarity/
├── README.md
├── benchmark_similarity_status.csv
├── benchmark_similarity_summary.csv
├── density/
│   ├── density_group_definitions.csv
│   ├── fixed_similarity_thresholds.csv
│   ├── model_specific_quantile_boundaries.csv
│   └── success_similarity_density_curves.csv
└── results/
    ├── boltz1/
    │   └── <PDB_ID>/
    └── boltz2/
        └── <PDB_ID>/
```

## Top-level tables

### `benchmark_similarity_status.csv`

Records the similarity-analysis status of all 100 benchmark complexes for
each temporal reference universe. It identifies complexes with per-target reports and complexes for which no eligible final PLINDER system was available.

### `benchmark_similarity_summary.csv`

Contains the query-level data used to analyse structural prediction success
as a function of similarity and training-set density.

Each record associates a benchmark complex, model and similarity metric
with fields such as:

- prediction success;
- highest observed similarity;
- number of neighbours above the selected threshold;
- mean neighbour similarity;
- assigned density group;
- similarity quantile;
- metric-specific threshold.

## Per-complex results

Each report-bearing directory under `results/boltz1/` or `results/boltz2/`
contains five files:

### `pairwise_scores.parquet`

Pairwise similarity values between the benchmark query and candidate
PLINDER reference systems.

### `pdb_target_scores.parquet`

Similarity values aggregated at the target PDB-entry level.

### `similarity_summary.csv`

Compact summary of the six similarity dimensions, including:

- number of reference systems considered;
- number and percentage of systems successfully evaluated;
- minimum, median and maximum calculated similarity;
- number of systems above the selected threshold;
- percentage of systems above the threshold;
- minimum, median and maximum passing similarity.

### `similarity_summary.parquet`

Parquet representation of the same summary information.

### `similarity_report.txt`

Human-readable report for manual inspection.

## Density-analysis tables

### `density/density_group_definitions.csv`

Definitions and observed ranges of the no-neighbour, low-density and
high-density training-reference groups.

### `density/fixed_similarity_thresholds.csv`

Metric-specific similarity thresholds used to define neighbouring
reference systems.

These thresholds were selected in a preceding analysis and then kept fixed
for the 100-complex benchmark. They were not re-optimised using prediction
success on this benchmark.

### `density/model_specific_quantile_boundaries.csv`

Similarity-quantile boundaries calculated separately for the Boltz-1 and
Boltz-2 reference universes.

### `density/success_similarity_density_curves.csv`

Aggregated prediction-success values used to generate the
success-versus-similarity and training-set-density figures.

## Data availability

The repository includes compact summaries and detailed per-query result
tables. It does not distribute:

- downloaded raw mmCIF structures;
- materialised PLINDER systems;
- Foldseek or MMseqs2 databases;
- container images;
- raw Boltz prediction directories.

This keeps the repository manageable while retaining the processed data
needed to inspect and reproduce the reported similarity analyses.

## Reproducibility

Scripts used for preparing and consolidating these results are documented
under [`../scripts/similarity/`](../scripts/similarity/).

The analysis notebook reads these processed files and reproduces the
reported association between prediction success, maximum training-reference
similarity and local training-reference density.
