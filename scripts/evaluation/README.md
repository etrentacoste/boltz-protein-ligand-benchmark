# Structural evaluation with OpenStructure

This directory contains the scripts used to evaluate Boltz-1 and Boltz-2
protein–ligand predictions against the experimental benchmark structures.

For each of the five predicted poses per model, the workflow:

1. creates an explicit ligand SDF from the Boltz prediction;
2. removes water and non-selected non-polymer components from the experimental
   reference structure;
3. calculates ligand RMSD and lDDT-PLI with OpenStructure;
4. summarises pose-level results into one model-level result per complex.

The scripts require OpenStructure and Gemmi. The study used OpenStructure
2.11.1.

## Scripts

- `evaluate_complex.py`: evaluates the five Boltz-1 and five Boltz-2 poses for
  one benchmark complex.
- `make_boltz_ligand_sdf.py`: creates explicit ligand SDF files from Boltz NPZ
  topology and predicted mmCIF coordinates.
- `summarize_boltz_results.py`: converts pose-level OpenStructure output into
  per-complex result tables.

## Inputs

The runtime directory passed through `--root` must contain one directory per
PDB entry. For a PDB identifier `<PDB_ID>`, the expected files include:

```text
<root>/<PDB_ID>/experimental_<PDB_ID>.cif
<root>/<PDB_ID>/production/boltz1/<PDB_ID>_model_0.cif
...
<root>/<PDB_ID>/production/boltz1/<PDB_ID>_model_4.cif
<root>/<PDB_ID>/production/boltz2/<PDB_ID>_model_0.cif
...
<root>/<PDB_ID>/production/boltz2/<PDB_ID>_model_4.cif
```

The published benchmark configurations are in
[`benchmark/configs/`](../../benchmark/configs/), with one file named
`<PDB_ID>.json` per complex.

## Run one complex

```bash
python scripts/evaluation/evaluate_complex.py <PDB_ID> \
    --root /path/to/runtime_complex_directory \
    --config-dir benchmark/configs
```

For example:

```bash
python scripts/evaluation/evaluate_complex.py 9VCZ \
    --root /path/to/boltz_visualization \
    --config-dir benchmark/configs
```

If `--config-dir` is omitted, the script instead reads:

```text
<root>/<PDB_ID>/config.json
```

This preserves compatibility with the original local runtime layout.

## Success definition

A pose was considered successful only when both conditions were satisfied:

```text
BiSyRMSD <= 2.0 Å
lDDT-PLI >= 0.8
```

For each model and complex, `top-1` refers to the model-ranked first pose and
`best-of-five` refers to whether at least one of the five generated poses met
both thresholds.

Experimental structures were used only for retrospective evaluation; they were
not used to choose the model-ranked top-1 prediction.
