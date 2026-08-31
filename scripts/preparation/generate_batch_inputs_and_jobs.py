#!/usr/bin/env python3

import argparse
import csv
import json
from pathlib import Path

import yaml


def validate_msa(path, expected_sequence):
    if not path.is_file():
        raise FileNotFoundError(path)

    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise RuntimeError(f"Empty MSA: {path}")

    query = rows[0]["sequence"].replace("-", "").upper()

    if query != expected_sequence.upper():
        raise RuntimeError(
            f"MSA query mismatch: {path}"
        )


def build_yaml(config, model):
    sequences = []

    for protein in config["proteins"]:
        msa_path = Path(protein["msa"])

        validate_msa(
            msa_path,
            protein["sequence"],
        )

        sequences.append({
            "protein": {
                "id": protein["id"],
                "sequence": protein["sequence"],
                "msa": str(msa_path),
            }
        })

    for ligand in config["ligands"]:
        ligand_entry = {
            "id": ligand["id"],
        }

        if ligand.get("smiles"):
            ligand_entry["smiles"] = ligand["smiles"]
        elif ligand.get("ccd"):
            ligand_entry["ccd"] = ligand["ccd"]
        else:
            raise RuntimeError(
                "Ligand has neither SMILES nor CCD"
            )

        sequences.append({
            "ligand": ligand_entry,
        })

    for cofactor in config.get("cofactors", []):
        if isinstance(cofactor, str):
            cofactor = {
                "id": f"C{len(sequences)}",
                "ccd": cofactor,
            }

        entry = {"id": cofactor["id"]}

        if cofactor.get("smiles"):
            entry["smiles"] = cofactor["smiles"]
        else:
            entry["ccd"] = cofactor["ccd"]

        sequences.append({"ligand": entry})

    for metal in config.get("metals", []):
        if isinstance(metal, str):
            metal = {
                "id": f"M{len(sequences)}",
                "ccd": metal,
            }

        sequences.append({
            "ligand": {
                "id": metal["id"],
                "ccd": metal["ccd"],
            }
        })

    result = {
        "version": 1,
        "sequences": sequences,
    }

    if model == "boltz2":
        binders = [
            ligand["id"]
            for ligand in config["ligands"]
            if ligand.get("affinity_binder", False)
        ]

        if len(binders) != 1:
            raise RuntimeError(
                "Boltz-2 requires exactly one "
                "affinity binder"
            )

        result["properties"] = [
            {
                "affinity": {
                    "binder": binders[0],
                }
            }
        ]

    return result


def total_protein_length(config):
    total = 0

    for protein in config["proteins"]:
        copies = (
            len(protein["id"])
            if isinstance(protein["id"], list)
            else 1
        )

        total += len(protein["sequence"]) * copies

    return total


def resources(config):
    length = total_protein_length(config)

    if length <= 400:
        memory = "48G"
        pilot_time = "02:00:00"
        production_time = "04:00:00"
    elif length <= 800:
        memory = "64G"
        pilot_time = "03:00:00"
        production_time = "06:00:00"
    else:
        memory = "96G"
        pilot_time = "04:00:00"
        production_time = "08:00:00"

    return {
        "cpus": 8,
        "memory": memory,
        "pilot_time": pilot_time,
        "production_time": production_time,
        "protein_length": length,
    }


def slurm_header(
    job_name,
    memory,
    time_limit,
    stdout,
    stderr,
    root,
):
    return f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem={memory}
#SBATCH --time={time_limit}
#SBATCH --chdir={root}
#SBATCH --output={stdout}
#SBATCH --error={stderr}

set -euo pipefail

module purge
module load conda/python3
eval "$(conda shell.bash hook)"
conda activate "${{BOLTZ_CONDA_ENV:?Set BOLTZ_CONDA_ENV to the Boltz conda environment path}}"

echo "Node: $(hostname)"
echo "Start: $(date)"
echo "Job ID: $SLURM_JOB_ID"
which python
which boltz
python --version

"""


def pilot_script(
    pdb_id,
    yaml_path,
    output_dir,
    resource,
    root,
):
    logs = root / "logs" / "pilot"
    logs.mkdir(parents=True, exist_ok=True)

    header = slurm_header(
        f"pilot_b2_{pdb_id}",
        resource["memory"],
        resource["pilot_time"],
        logs / f"boltz2_{pdb_id}_%j.out",
        logs / f"boltz2_{pdb_id}_%j.err",
        root,
    )

    command = f"""boltz predict {yaml_path} \\
    --model boltz2 \\
    --accelerator cpu \\
    --devices 1 \\
    --recycling_steps 1 \\
    --sampling_steps 20 \\
    --diffusion_samples 1 \\
    --max_parallel_samples 1 \\
    --sampling_steps_affinity 20 \\
    --diffusion_samples_affinity 1 \\
    --seed 12345 \\
    --override \\
    --out_dir {output_dir}

echo "Finish: $(date)"
"""

    return header + command


def production_script(
    pdb_id,
    model,
    yaml_path,
    output_dir,
    resource,
    root,
):
    logs = root / "logs" / "production"
    logs.mkdir(parents=True, exist_ok=True)

    short_model = "b1" if model == "boltz1" else "b2"

    header = slurm_header(
        f"prod_{short_model}_{pdb_id}",
        resource["memory"],
        resource["production_time"],
        logs / f"{short_model}_{pdb_id}_%j.out",
        logs / f"{short_model}_{pdb_id}_%j.err",
        root,
    )

    affinity_options = ""

    if model == "boltz2":
        affinity_options = """ \\
    --sampling_steps_affinity 200 \\
    --diffusion_samples_affinity 5"""

    command = f"""boltz predict {yaml_path} \\
    --model {model} \\
    --accelerator cpu \\
    --devices 1 \\
    --recycling_steps 3 \\
    --sampling_steps 200 \\
    --diffusion_samples 5 \\
    --max_parallel_samples 1{affinity_options} \\
    --seed 12345 \\
    --override \\
    --out_dir {output_dir}

echo "Finish: $(date)"
"""

    return header + command


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--batch",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--root",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--report",
        type=Path,
        required=True,
    )

    args = parser.parse_args()

    with args.batch.open(newline="") as handle:
        batch_rows = list(csv.DictReader(handle))

    generated_dir = (
        args.root
        / "scripts"
        / "generated"
    )
    generated_dir.mkdir(parents=True, exist_ok=True)

    results = []

    for row in batch_rows:
        pdb_id = row["pdb_id"].upper()
        config_path = (
            args.root
            / "config"
            / f"{pdb_id}.json"
        )
        config = json.loads(config_path.read_text())
        resource = resources(config)

        yaml_paths = {}

        for model in ["boltz1", "boltz2"]:
            yaml_data = build_yaml(config, model)
            yaml_path = (
                args.root
                / "inputs"
                / model
                / f"{pdb_id}.yaml"
            )
            yaml_path.parent.mkdir(
                parents=True,
                exist_ok=True,
            )
            yaml_path.write_text(
                yaml.safe_dump(
                    yaml_data,
                    sort_keys=False,
                )
            )
            yaml_paths[model] = yaml_path

        pilot_path = (
            generated_dir
            / f"pilot_{pdb_id}_boltz2.slurm"
        )
        pilot_output = (
            args.root
            / "results"
            / "pilot"
            / "boltz2"
            / pdb_id
        )

        pilot_path.write_text(
            pilot_script(
                pdb_id,
                yaml_paths["boltz2"],
                pilot_output,
                resource,
                args.root,
            )
        )
        pilot_path.chmod(0o755)

        production_paths = {}

        for model in ["boltz1", "boltz2"]:
            production_path = (
                generated_dir
                / f"production_{pdb_id}_{model}.slurm"
            )
            production_output = (
                args.root
                / "results"
                / "production"
                / model
                / pdb_id
            )

            production_path.write_text(
                production_script(
                    pdb_id,
                    model,
                    yaml_paths[model],
                    production_output,
                    resource,
                    args.root,
                )
            )
            production_path.chmod(0o755)
            production_paths[model] = production_path

        results.append({
            "pdb_id": pdb_id,
            "protein_length": (
                resource["protein_length"]
            ),
            "memory": resource["memory"],
            "pilot_time": (
                resource["pilot_time"]
            ),
            "production_time": (
                resource["production_time"]
            ),
            "boltz1_yaml": str(
                yaml_paths["boltz1"]
            ),
            "boltz2_yaml": str(
                yaml_paths["boltz2"]
            ),
            "pilot_script": str(pilot_path),
            "boltz1_script": str(
                production_paths["boltz1"]
            ),
            "boltz2_script": str(
                production_paths["boltz2"]
            ),
            "status": "SUCCESS",
        })

        print(
            pdb_id,
            f"length={resource['protein_length']}",
            f"memory={resource['memory']}",
            "SUCCESS",
        )

    with args.report.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=results[0].keys(),
        )
        writer.writeheader()
        writer.writerows(results)

    print()
    print("Complexes:", len(results))
    print("Generated:", len(results))
    print("Report:", args.report)


if __name__ == "__main__":
    main()
