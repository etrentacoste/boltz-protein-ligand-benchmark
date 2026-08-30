# Boltz-2 inputs

This directory contains representative YAML inputs used for Boltz-2
protein–ligand structure and affinity prediction.

The example `31WA.yaml` specifies:

- the curated protein sequence;
- the protein-chain identifier;
- the precomputed MSA path;
- the selected ligand as a SMILES string;
- the ligand identifier used as the Boltz-2 affinity binder.

The absolute MSA path records the location used during the original
calculations on Iridis. Users reproducing the workflow must replace this
path with the location of their own MSA file.

MSA files are not stored in this GitHub repository because of their size.
