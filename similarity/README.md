# PLINDER similarity analysis results

This directory contains processed tables generated during the analysis of
prediction success as a function of similarity to the PLINDER training set.

The analysis considered six similarity dimensions:

- protein sequence similarity;
- protein structure similarity;
- binding-pocket similarity;
- protein–ligand interaction similarity;
- ligand-pose similarity;
- ligand-chemistry similarity.

PLINDER release `2024-06/v2` was used. For the common comparison between
prediction methods, training systems were restricted using the conservative
cutoff date of 30 September 2021.

## Directory structure

### `processed/`

`proper_ligand_similarity_and_prediction_results.csv` contains the combined
prediction and similarity results after ligand-level filtering.

`common_complete_case_similarity_and_prediction_results.csv` contains the
complete-case dataset used for comparisons across prediction methods and
similarity dimensions.

### `thresholds/`

`candidate_threshold_effects.csv` reports the effects obtained for candidate
neighbour thresholds.

`selected_similarity_thresholds.csv` records the threshold selected for each
similarity metric.

`selected_threshold_bootstrap_summary.csv` contains bootstrap estimates for
the selected thresholds.

The bootstrap procedure used 2,000 query-level replicates with random seed 42.

### `density/`

Files ending in `_density_query_table.csv` contain the query-level observations
used in the density analyses.

Files ending in `_density_curve_summary.csv` contain the aggregated values
plotted in the success-versus-similarity figures.

Training-set density was divided into three groups:

- no neighbours;
- low density;
- high density.

The numerical boundaries of these groups depend on the selected similarity
metric and are recorded in the corresponding output tables.

## Large source files

The following source files are not stored in this GitHub repository because
of their size:

- `all_similarity_scores.parquet`;
- `predictions/af3.csv`;
- `predictions/boltz.csv`;
- `predictions/chai.csv`;
- `predictions/protenix.csv`.

These files were obtained or derived from the PLINDER benchmark data and are
used by the preparation notebook. The smaller processed tables required to
inspect the reported analyses are provided in this directory.

The analysis workflow is documented in:

- `notebooks/01_plinder_similarity_preparation.ipynb`;
- `notebooks/02_plinder_similarity_analysis.ipynb`.
