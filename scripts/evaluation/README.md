# Structural evaluation

This directory contains the scripts used to evaluate Boltz-1 and Boltz-2
protein–ligand complex predictions against experimentally determined
reference structures.

## Available scripts

### `evaluate_complex.py`

Main evaluation workflow for one benchmark system.

The script:

1. reads the curated system configuration from `config.json`;
2. loads the experimental mmCIF reference structure;
3. retains the selected protein chains;
4. removes water molecules and unrelated non-polymer components;
5. retains the selected ligand, cofactors and metals;
6. extracts the ligand from each predicted Boltz structure;
7. evaluates five Boltz-1 and five Boltz-2 structural samples;
8. calls OpenStructure `compare-ligand-structures`;
9. calculates symmetry-aware ligand RMSD;
10. calculates lDDT-PLI;
11. records ligand coverage and binding-site mapping information;
12. validates that a ligand assignment was produced;
13. calls `summarize_boltz_results.py` to classify and summarise the results.

Example:

```bash
conda activate ost-env

python scripts/evaluation/evaluate_complex.py 9VCZ \
    --root /path/to/boltz_visualization
