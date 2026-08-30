# Boltz protein–ligand benchmark

Code, processed results and analysis workflows used to benchmark Boltz-1
and Boltz-2 for protein–ligand complex structure prediction.

## Overview

This repository accompanies a research project evaluating the performance
of Boltz-1 and Boltz-2 on recently released experimentally determined
protein–ligand complexes.

The project includes:

- construction and manual curation of a 100-complex benchmark;
- preparation of protein, ligand, MSA and cofactor inputs;
- protein–ligand structure prediction with Boltz-1 and Boltz-2;
- evaluation of five structural samples per model and complex;
- calculation of ligand RMSD and lDDT-PLI using OpenStructure;
- classification of top-1 and best-of-five prediction success;
- comparison of Boltz-1 and Boltz-2 performance;
- PLINDER-based analysis of similarity to training-set structures;
- analysis of training-set neighbourhood density;
- bootstrap estimation of uncertainty;
- molecular visualisation with PyMOL.

The repository is being organised alongside the final project report.
Additional scripts, processed tables and documentation will be added
progressively.

## Benchmark composition

The final benchmark contains 100 experimentally determined
protein–ligand complexes.

It combines:

- 30 systems selected during the initial manual curation stage;
- 70 additional systems selected using an automated RCSB screening,
  structural-quality filtering, ligand-validation filtering and
  target-diversity workflow.

The selection procedure considered:

- experimental method and structural resolution;
- crystallographic R-free;
- ligand occupancy, RSCC and RSR;
- ligand chemical composition;
- covalent bonding;
- relevant protein chains;
- repeated protein chains;
- multiple ligand instances;
- cofactors and metals;
- crystallographic solvents and artefacts;
- diversity of protein targets.

Experimental structures are identified by their RCSB PDB accession codes.
The original RCSB coordinate files are not redistributed in this
repository.

## Prediction protocol

Each benchmark system was predicted independently with Boltz-1 and
Boltz-2.

The production inference settings were:

| Parameter | Value |
|---|---:|
| Recycling steps | 3 |
| Sampling steps | 200 |
| Diffusion samples | 5 |
| Maximum parallel samples | 1 |
| Random seed | 12345 |
| Boltz-2 affinity sampling steps | 200 |
| Boltz-2 affinity diffusion samples | 5 |
| Accelerator used in this study | CPU |
| Devices per job | 1 |

Five structural samples were generated for each model and benchmark
system. Boltz-2 was additionally used to generate affinity outputs.

Calculations were executed as Slurm jobs on the Iridis high-performance
computing cluster at the University of Southampton.

Example inference scripts and input YAML files are provided under
`scripts/inference/` and `inputs/`.

## Structural evaluation

Predicted complexes were evaluated against their experimentally
determined structures using OpenStructure.

The evaluation workflow:

1. retained the curated protein chains and selected experimental ligand;
2. removed water molecules and unrelated crystallographic components;
3. retained explicitly selected cofactors and metals;
4. extracted the ligand from every predicted complex;
5. mapped the predicted and experimental protein chains;
6. superposed the binding sites;
7. calculated symmetry-aware ligand RMSD;
8. calculated lDDT-PLI;
9. recorded ligand coverage and binding-site mapping information;
10. classified top-1 and best-of-five success.

A pose was classified as successful when both of the following conditions
were satisfied:

- ligand RMSD < 2 Å;
- lDDT-PLI > 0.8.

Failure to obtain a valid ligand assignment was treated as an evaluation
error and was not automatically classified as an unsuccessful pose.

## Benchmark results

Across the 100 benchmark complexes, the aggregate results were:

| Model | Top-1 success | Best-of-five success |
|---|---:|---:|
| Boltz-1 | 39/100 (39%) | 43/100 (43%) |
| Boltz-2 | 42/100 (42%) | 46/100 (46%) |

Detailed pose-level and model-level results are provided under
`results/`.

## PLINDER similarity analysis

Similarity to experimentally determined complexes in PLINDER was
evaluated using multiple complementary representations:

- protein sequence similarity;
- protein structure similarity;
- binding-pocket similarity;
- protein–ligand interaction similarity;
- ligand-pose similarity;
- ligand-chemistry similarity.

The analysis used the PLINDER `2024-06/v2` data release.

For each query system and similarity representation, the analysis
identified related training-set systems above a predefined similarity
threshold. Systems were then divided into:

- no training-set neighbours;
- low training-set density;
- high training-set density.

The low- and high-density groups were separated using the median positive
neighbour count for each similarity representation.

Uncertainty was estimated using 2,000 query-level bootstrap replicates
with random seed 42.

The similarity-processing and statistical-analysis code will be provided
under `scripts/similarity/` and `notebooks/`.

## Repository structure

```text
benchmark/     Benchmark composition and selection tables
configs/       Curated system configurations
inputs/        Boltz-1 and Boltz-2 YAML inputs
scripts/       Selection, preparation, inference and evaluation code
results/       Processed pose-level and model-level results
notebooks/     Final analysis notebooks
figures/       Figures used in the report
environments/  Software environment specifications
docs/          Extended methodological documentation
