# Boltz-1 inputs

This directory contains representative YAML inputs used for Boltz-1
structure prediction.

The example `31WA.yaml` specifies:

- the curated protein sequence;
- the protein-chain identifier;
- the precomputed MSA path;
- the selected ligand as a SMILES string.

The absolute MSA path records the location used during the original
calculations on Iridis. Users reproducing the workflow must replace this
path with the location of their own MSA file.

MSA files are not stored in this GitHub repository because of their size.
