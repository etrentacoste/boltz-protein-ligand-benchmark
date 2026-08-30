# Processed benchmark results

This directory contains processed outputs from the Boltz-1 and Boltz-2
benchmark.

Available files include:

- `benchmark_100_model_summary.csv`: one row per model and complex;
- `benchmark_100_pose_results.csv`: one row per predicted pose;
- `benchmark_100_affinity_results.csv`: processed Boltz-2 affinity outputs;
- `benchmark_performance_summary.csv`: aggregate top-1 and best-of-five
  success rates;
- `similarity_density_summary.csv`: processed PLINDER similarity and
  neighbourhood-density results;
- `bootstrap_summary.csv`: bootstrap estimates and confidence intervals.

The primary structural success criterion was:

- ligand RMSD < 2 Å; and
- lDDT-PLI > 0.8.

Aggregate benchmark performance:

| Model | Top-1 success | Best-of-five success |
|---|---:|---:|
| Boltz-1 | 39/100 (39%) | 43/100 (43%) |
| Boltz-2 | 42/100 (42%) | 46/100 (46%) |

Large raw prediction files are not stored in this GitHub repository.
