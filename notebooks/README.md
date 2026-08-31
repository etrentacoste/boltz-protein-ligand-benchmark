# Benchmark analysis notebook

## `boltz_benchmark_analysis.ipynb`

This notebook contains the analysis developed for the 100-complex
Boltz-1 and Boltz-2 protein–ligand benchmark.

The notebook:

- loads the structural prediction results for Boltz-1 and Boltz-2;
- applies the predefined pose-success criterion;
- combines prediction outcomes with PLINDER similarity measurements;
- applies similarity thresholds established independently from the benchmark;
- calculates local training-set density;
- divides benchmark systems into similarity quantiles;
- compares Boltz-1 and Boltz-2 performance;
- generates the numerical results and figures reported in the dissertation.

A predicted pose was classified as successful when:

- BiSyRMSD < 2 Å; and
- lDDT-PLI > 0.8.

## Inputs

The notebook requires:

- the 100-complex benchmark results;
- pairwise PLINDER similarity scores;
- the previously selected similarity thresholds.

Large raw prediction and PLINDER files are not stored in this repository.
Processed tables required to inspect the reported results are available
under `results/` and `similarity/`.

## Data location

Set the `BOLTZ_BENCHMARK_DATA` environment variable to the directory
containing the external input files:

```bash
export BOLTZ_BENCHMARK_DATA=/path/to/analysis/data
```
