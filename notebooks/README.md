# Analysis notebooks

This directory contains the Jupyter notebooks used to analyse the
relationship between Boltz prediction success and similarity to the
PLINDER training set.

## Notebooks

### `01_plinder_similarity_preparation.ipynb`

Prepares and validates the combined benchmark, prediction-success and
PLINDER-similarity data used in the downstream analysis.

The notebook includes data loading, column validation, similarity metric
preparation and construction of analysis-ready tables.

### `02_plinder_similarity_analysis.ipynb`

Performs the final similarity and training-set-density analyses and
generates the corresponding figures.

The analysed similarity dimensions include:

- protein sequence similarity;
- protein structure similarity;
- binding-pocket similarity;
- protein–ligand interaction similarity;
- ligand-pose similarity;
- ligand-chemistry similarity.

Training-set density was defined from the number of PLINDER neighbours
above a metric-specific similarity threshold. Systems were grouped as
having no neighbours, low training-set density or high training-set
density.

Uncertainty was estimated using 2,000 query-level bootstrap replicates
with random seed 42.

## Execution order

Run the notebooks in numerical order:

1. `01_plinder_similarity_preparation.ipynb`
2. `02_plinder_similarity_analysis.ipynb`

The notebooks are provided with their stored outputs to preserve the
figures and numerical summaries used in the project report.

## Data availability

Large PLINDER source files and cache directories are not included in
this repository. The notebooks therefore require a local PLINDER
installation and access to the relevant PLINDER release.

The analysis used PLINDER release `2024-06/v2`.
