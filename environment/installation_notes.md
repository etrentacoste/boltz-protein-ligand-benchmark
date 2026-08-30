# Installation notes

The commands below document the main packages required by the workflow. They
are not intended to recreate the University of Southampton HPC configuration
or Slurm installation.

Separate environments are recommended for inference, structural evaluation,
and similarity analysis.

## Boltz inference environment

```bash
conda create -n boltz-inference python=3.11
conda activate boltz-inference
python -m pip install "boltz==2.2.1"
