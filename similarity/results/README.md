# Per-complex similarity results

This directory contains the similarity results between each eligible benchmark
complex and the PLINDER-based reference universe defined for each Boltz model.

```text
results/
├── boltz1/<PDB_ID>/
└── boltz2/<PDB_ID>/
```

Each model directory contains one subdirectory for every benchmark complex for
which a final eligible PLINDER query system was available.

There are 90 complex directories for Boltz-1 and 90 for Boltz-2. The remaining
10 benchmark complexes did not yield an eligible final PLINDER query system;
their status is recorded in
[`../benchmark_similarity_status.csv`](../benchmark_similarity_status.csv).

## Files for each complex

Each `<PDB_ID>/` directory contains:

- `pairwise_scores.parquet`: similarity values between the benchmark query
  system and individual PLINDER reference systems;
- `pdb_target_scores.parquet`: pairwise values aggregated by target PDB entry;
- `similarity_summary.csv`: compact summary of the six similarity metrics;
- `similarity_summary.parquet`: the same summary in Parquet format;
- `similarity_report.txt`: human-readable report for the complex.

The six summary metrics are protein sequence similarity, protein structure
similarity, pocket similarity, interaction similarity, ligand-pose similarity
and ligand-chemistry similarity.

The calculation workflow is documented in
[`../../scripts/similarity/`](../../scripts/similarity/). Dataset-level
summaries and the density-analysis tables are located one directory above in
[`../`](../).
