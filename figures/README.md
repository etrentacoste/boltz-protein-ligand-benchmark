# Figures

This directory contains selected final figures supporting the benchmark and
similarity analyses. Preliminary versions, duplicate renderings, individual
PyMOL sessions, and intermediate panels are not included.

## `benchmark/`

`figure_model_comparison.png` summarises the benchmark-level comparison
between Boltz-1 and Boltz-2.

`figure_failure_mechanisms.png` summarises the principal categories of
prediction outcome and failure observed during structural evaluation.

## `similarity/`

The six `*_boltz1_boltz2_density_lines_bootstrap_ci.pdf` files show prediction
success across similarity quantiles and training-set-density groups for:

- protein sequence;
- protein structure;
- binding pocket;
- protein–ligand interactions;
- ligand pose;
- ligand chemistry.

Points show mean prediction success. Shaded regions represent pointwise 95%
confidence intervals obtained using 2,000 query-level bootstrap replicates
with random seed 42. Sample-size labels report the number of benchmark
complexes contributing to each point.

`candidate_threshold_bootstrap_thesis_A4.pdf` shows the bootstrap analysis used
to select the neighbour threshold for each similarity metric.

The underlying numerical tables are available in
[`../similarity/`](../similarity/).

## `selection/`

The three PyMOL figures illustrate examples from benchmark curation:

- `figure_1_include_9RQV_v2.png`: a suitable protein–ligand complex retained
  for the benchmark;
- `figure_2_manual_review_24AP_v2.png`: a system requiring manual review
  because of its metal-, nucleotide-, and oligomer-dependent context;
- `figure_3_exclude_27FM_v3.png`: a system excluded because its small
  components are metal ions and crystallographic additives rather than a
  conventional organic ligand.

These figures are illustrative examples of the curation logic rather than a
complete representation of every filtering rule.

## `concepts/`

`9RQV_docking_concept.png` illustrates conventional docking as placement of a
ligand into a predefined binding site on a supplied receptor structure.

`9RQV_cofolding_samples_v3.png` illustrates protein–ligand cofolding using five
complete structural samples generated for the same system. In cofolding, the
protein and ligand coordinates are predicted jointly rather than placing the
ligand into a rigid, predefined receptor structure.

## Colours

Unless otherwise stated in an individual figure:

- grey represents the protein structure;
- cyan highlights binding-pocket residues or surfaces;
- distinct ligand colours distinguish experimental and predicted poses.
