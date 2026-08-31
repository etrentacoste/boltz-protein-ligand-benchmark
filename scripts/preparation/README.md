# Simulation-input preparation

This directory contains the original scripts used to prepare the 100 benchmark
complexes for Boltz-1 and Boltz-2 inference.

The scripts prepare the per-complex runtime layout, extract the required
protein and ligand information from the benchmark configuration, generate or
reuse MSAs, write Boltz YAML inputs, retain selected cofactors and metals, and
create Slurm job scripts for batch execution.

The raw structures, generated MSAs, production predictions and Slurm logs are
not distributed in this repository. They can be recreated from the benchmark
tables, configurations and scripts provided here.

## Scripts

| Script | Role |
|---|---|
| `prepare_simulation_batch.py` | Prepares the runtime files and per-complex configuration for one simulation batch. |
| `generate_batch_msas.py` | Generates or reuses MSAs for all complexes in one batch. |
| `generate_batch_inputs_and_jobs.py` | Generates Boltz input YAML files and batch-specific Slurm jobs. |
| `generate_boltz_inputs.py` | Generates a Boltz input YAML file from one complex configuration. |
| `generate_slurm_jobs.py` | Generates Slurm jobs for one complex configuration. |
| `generate_msa_only.sh` | Shell helper for MSA generation. |
| `generate_msa_with_retries.sh` | Shell helper that retries failed MSA-generation attempts. |

## Typical batch workflow

The batch-oriented scripts require a runtime root containing the benchmark
complex directories and a batch number:

```bash
python scripts/preparation/prepare_simulation_batch.py \
    --batch <BATCH_NUMBER> \
    --root /path/to/runtime_root \
    --report /path/to/preparation_report.csv

python scripts/preparation/generate_batch_msas.py \
    --batch <BATCH_NUMBER> \
    --root /path/to/runtime_root \
    --report /path/to/msa_report.csv

python scripts/preparation/generate_batch_inputs_and_jobs.py \
    --batch <BATCH_NUMBER> \
    --root /path/to/runtime_root \
    --report /path/to/input_and_job_report.csv
```

The generated YAML files can be inspected against the model-specific examples
in [`../../inputs/`](../../inputs/). The generated Slurm jobs run the Boltz
commands documented in [`../inference/`](../inference/).

## Slurm configuration

`generate_batch_inputs_and_jobs.py` writes a Slurm template that requires the
environment variable `BOLTZ_CONDA_ENV`. Define and export it before submitting
the generated jobs:

```bash
export BOLTZ_CONDA_ENV="/absolute/path/to/boltz-conda-environment"
```

The supplied template was used on the Iridis cluster and includes its
`module load conda/python3` command and `batch` partition. Adapt module names,
partition, resources and scheduler directives to the target system before
running it elsewhere.

## Cohort configuration

The original batch workflow used the published complex configurations in
[`../../benchmark/configs/`](../../benchmark/configs/). Each configuration
specifies the selected protein chains, ligand, cofactors and metals retained
when preparing the model input and experimental evaluation reference.
