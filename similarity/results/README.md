# Prediction and evaluation results

This directory contains the consolidated numerical results obtained by
predicting the 100 benchmark protein–ligand complexes with Boltz-1 and
Boltz-2 and evaluating the predicted poses against their experimental
structures.

The structural predictions themselves are not included in this repository
because the complete set of predicted mmCIF files is substantially larger
than the processed result tables. The commands and scripts required to
reproduce the predictions and evaluations are provided under
[`scripts/inference/`](../scripts/inference/) and
[`scripts/evaluation/`](../scripts/evaluation/).

## Files

### `benchmark_100_pose_results.csv`

Pose-level evaluation table.

Each of the 100 benchmark complexes was predicted with Boltz-1 and Boltz-2.
Five diffusion samples were generated per model, producing 1,000 evaluated
poses in total:

- 100 benchmark complexes;
- 2 prediction models;
- 5 structural samples per model.

The table includes the ligand BiSyRMSD, lDDT-PLI, ligand coverage and other
OpenStructure-derived quantities for each individual pose.

### `benchmark_100_model_summary.csv`

Model-level summary containing one row for every benchmark complex and
prediction model, for a total of 200 rows.

For each complex and model, the table reports:

- top-1 success;
- best-of-five success;
- number of successful poses among the five samples;
- top-1 BiSyRMSD;
- top-1 lDDT-PLI;
- best BiSyRMSD among the five samples;
- best lDDT-PLI among the five samples;
- Boltz-2 affinity outputs when available.

The top-1 structure corresponds to the model-ranked prediction rather than
the pose selected retrospectively using the experimental structure.

### `benchmark_100_affinity_predictions.csv`

Boltz-2 affinity outputs for the 100 benchmark complexes.

The file reports:

- `affinity_pred_value`: the affinity-related scalar produced by Boltz-2;
- `binder_probability`: the Boltz-2 predicted probability that the supplied
  ligand is a binder.

These quantities were retained as model outputs and were not used to define
structural prediction success.

## Structural evaluation

Predicted complexes were compared with the experimentally determined
structures using OpenStructure 2.11.1.

The principal structural metrics were:

- **BiSyRMSD**, measuring ligand pose deviation after binding-site-aware
  structural superposition;
- **lDDT-PLI**, measuring preservation of local protein–ligand contacts;
- **ligand coverage**, measuring the fraction of the ligand included in the
  structural assignment.

A pose was classified as successful only when both conditions were met:

```text
BiSyRMSD <= 2.0 Å
lDDT-PLI >= 0.8
```
A complex was classified as:
- top-1 successful when its model-ranked first pose satisfied both
  thresholds;
- best-of-five successful when at least one of the five generated poses
  satisfied both thresholds.

Aggregate results
Across the 100-complex benchmark:
Model	Top-1 success	Best-of-five success
Boltz-1	39/100 (39%)	43/100 (43%)
Boltz-2	42/100 (42%)	46/100 (46%)


These values describe structural pose-prediction success on the benchmark
defined in [`benchmark/`](../benchmark/).
Reproducibility
The inference settings, example YAML inputs and execution commands are
documented under:
- [`scripts/inference/`](../scripts/inference/)
- [`inputs/`](../inputs/)
The OpenStructure evaluation and result-summarisation scripts are documented
under:
- [`scripts/evaluation/`](../scripts/evaluation/)
No experimental structure was used to select the model-ranked top-1
prediction. Experimental structures were used only during retrospective
evaluation.
