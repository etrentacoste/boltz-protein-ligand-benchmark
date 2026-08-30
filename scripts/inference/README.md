# Boltz inference

This directory contains representative Slurm scripts used to run Boltz-1
and Boltz-2 on the Iridis high-performance computing cluster.

## Available scripts

- `example_boltz1.slurm`: Boltz-1 production inference example.
- `example_boltz2.slurm`: Boltz-2 production structure and affinity
  inference example.

Both examples correspond to benchmark system `31WA`.

## Common inference settings

| Parameter | Value |
|---|---:|
| Accelerator | CPU |
| Devices | 1 |
| Recycling steps | 3 |
| Sampling steps | 200 |
| Diffusion samples | 5 |
| Maximum parallel samples | 1 |
| Random seed | 12345 |
| Override existing output | Enabled |

Boltz-2 additionally used:

| Parameter | Value |
|---|---:|
| Affinity sampling steps | 200 |
| Affinity diffusion samples | 5 |

Five diffusion samples were generated during a single model run. The five
poses were not generated using five independently varied random seeds.
The production run used one fixed seed, `12345`, to make inference
reproducible.

## Boltz-1 command

The essential command was:

```bash
boltz predict INPUT.yaml \
    --model boltz1 \
    --accelerator cpu \
    --devices 1 \
    --recycling_steps 3 \
    --sampling_steps 200 \
    --diffusion_samples 5 \
    --max_parallel_samples 1 \
    --seed 12345 \
    --override \
    --out_dir OUTPUT_DIRECTORY

## Boltz-2 command

The essential command was:

```boltz predict INPUT.yaml \
    --model boltz2 \
    --accelerator cpu \
    --devices 1 \
    --recycling_steps 3 \
    --sampling_steps 200 \
    --diffusion_samples 5 \
    --max_parallel_samples 1 \
    --sampling_steps_affinity 200 \
    --diffusion_samples_affinity 5 \
    --seed 12345 \
    --override \
    --out_dir OUTPUT_DIRECTORY

## Inputs and outputs

Representative YAML inputs are provided under:

`inputs/boltz1/`
`inputs/boltz2/`

Boltz-1 outputs include:

- five predicted complex structures in mmCIF format;
- one confidence JSON file per structural sample.

Boltz-2 outputs additionally include:

- affinity prediction JSON output;
- binder probability;
- affinity model scores.

## Paths

The Slurm examples retain the absolute paths used during the original
Iridis runs. These paths must be adapted when running the workflow in a
different environment.

## Computational environment

The calculations used:
- one Slurm task;
- eight CPU cores per task;
- CPU inference;
- memory allocated according to system size;
- the batch partition on Iridis.
Exact software versions will be documented under `environments/`.
