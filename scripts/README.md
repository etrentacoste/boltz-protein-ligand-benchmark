# Workflow scripts

This directory contains the scripts used to construct the benchmark, prepare
Boltz inputs, run inference, evaluate predictions and calculate similarity.

```text
selection/       RCSB querying and benchmark filtering
training_set/    Construction of PLINDER-based temporal reference manifests
preparation/     Complex preparation, ligand selection, MSA and YAML generation
inference/       Boltz-1 and Boltz-2 execution examples and Slurm submission
evaluation/      OpenStructure evaluation and success classification
similarity/      PLINDER, Foldseek, MMseqs2, SuCOS and density analysis
```

Each subdirectory has its own README with inputs, outputs and usage details.
