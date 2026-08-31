# Installation notes

The commands below document the principal software environments required by
the workflow. They are not intended to recreate the University of Southampton
Iridis HPC configuration, Apptainer installation or Slurm scheduler.

Separate environments are recommended for inference, similarity analysis,
structural evaluation and PyMOL visualisation.

## Boltz inference environment

```bash
conda create -n boltz-inference python=3.11
conda activate boltz-inference
python -m pip install "boltz==2.2.1"
```

## Similarity-analysis environment

```bash
conda create -n plinder-analysis python=3.10
conda activate plinder-analysis

conda install -c conda-forge \
    "duckdb=1.5.5" \
    "rdkit=2024.9.6" \
    pandas numpy scipy pyarrow matplotlib

conda install -c bioconda \
    "foldseek=10.941cd33" \
    "mmseqs2=18.8cc5c"

python -m pip install "plinder==0.2.25"
```

The SuCOS implementation used in this study is included directly in
[`../scripts/similarity/calculate_validate_sucos_shape.py`](../scripts/similarity/calculate_validate_sucos_shape.py).

## OpenStructure evaluation environment

```bash
mamba create -n ost-env \
    -c conda-forge \
    -c bioconda \
    python=3.12 \
    "openstructure=2.11.1" \
    gemmi
```

## PyMOL visualisation environment

```bash
mamba create -n pymol-env \
    -c conda-forge \
    "pymol-open-source=3.1.0"
```

## HPC components

The production similarity workflow used an Apptainer-compatible PLINDER
container and Slurm job submission on Iridis. These components are managed by
the HPC platform and should be loaded or configured according to the target
cluster documentation. Their exact platform versions were not recorded during
this study.
