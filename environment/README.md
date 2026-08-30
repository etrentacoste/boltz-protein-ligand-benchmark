# Software environments

This directory documents the principal software environments used for
inference, structural evaluation, and PLINDER similarity analysis.

Three separate Conda environments were used because Boltz, OpenStructure, and
PLINDER have different dependency requirements.

## Boltz inference

Boltz-1 and Boltz-2 predictions were generated with Boltz version `2.2.1`.

The two models were not installed as separate Python packages. They were
selected from the same installation using:

```text
--model boltz1
or
--model boltz2
```
Inference was performed on CPU nodes of the University of Southampton Iridis
HPC cluster. Jobs were submitted using Slurm.
The production calculations used:
- three recycling steps;
- 200 diffusion sampling steps;
- five diffusion samples;
- one sample processed at a time;
- fixed random seed 12345;
- one CPU device from the perspective of the Boltz command;
- eight allocated CPU cores per Slurm job.

Boltz-2 additionally used 200 affinity sampling steps and five affinity
diffusion samples.

## Structural evaluation
Structural evaluation was performed with OpenStructure version 2.11.1.
OpenStructure was used to calculate:
- binding-site-superposed ligand RMSD;
- lDDT-PLI;
- ligand coverage;
- binding-site backbone RMSD;
- lDDT-LP.

A prediction was classified as successful when both conditions were met:
ligand RMSD < 2 Å
lDDT-PLI > 0.8

## PLINDER similarity analysis
Similarity analyses used:
- PLINDER Python package version 0.2.25;
- PLINDER data release 2024-06/v2;
- 2,000 query-level bootstrap replicates;
- bootstrap random seed 42.

The complete PLINDER source data and pairwise similarity files are not
included in this repository because of their size. Processed analysis tables
are provided in [`../similarity/`](../similarity/).

Version record
Exact versions of the principal packages are recorded in
[`software_versions.txt`](software_versions.txt).
