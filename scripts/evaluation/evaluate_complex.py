import argparse
import json
import subprocess
import sys
from pathlib import Path

import gemmi


WATER_NAMES = {"HOH", "WAT", "DOD"}


def run(command):
    print()
    print("+", " ".join(str(value) for value in command))
    subprocess.run(
        [str(value) for value in command],
        check=True,
    )


def clean_reference(
    source,
    output,
    protein_chains,
    retained_nonpolymers,
):
    structure = gemmi.read_structure(str(source))
    model = structure[0]

    available_chains = [chain.name for chain in model]

    missing = [
        chain for chain in protein_chains
        if chain not in available_chains
    ]

    if missing:
        raise RuntimeError(
            f"Reference chains not found: {missing}. "
            f"Available chains: {available_chains}"
        )

    for chain_index in range(len(model) - 1, -1, -1):
        if model[chain_index].name not in protein_chains:
            del model[chain_index]

    ligand_counts = {}

    for chain in model:
        for residue_index in range(len(chain) - 1, -1, -1):
            residue = chain[residue_index]
            name = residue.name

            # Boltz predicts heavy atoms. Remove explicit hydrogen and
            # deuterium atoms from the experimental reference so that
            # ligand identity and graph matching use comparable entities.
            for atom_index in range(len(residue) - 1, -1, -1):
                atom = residue[atom_index]

                if atom.element.name in {"H", "D"}:
                    del residue[atom_index]

            if name in WATER_NAMES:
                del chain[residue_index]
                continue

            # Polymer residues are normally represented as ATOM records
            # (het_flag != H). Retain only explicitly selected HETATMs.
            if (
                residue.het_flag == "H"
                and name not in retained_nonpolymers
            ):
                del chain[residue_index]
                continue

            if name in retained_nonpolymers:
                ligand_counts[name] = (
                    ligand_counts.get(name, 0) + 1
                )

    output.parent.mkdir(parents=True, exist_ok=True)
    structure.name = output.stem
    structure.make_mmcif_document().write_file(str(output))

    print()
    print("Clean reference:", output)
    print("Retained chains:", [chain.name for chain in model])
    print("Retained nonpolymers:", ligand_counts)

    for required in retained_nonpolymers:
        if ligand_counts.get(required, 0) == 0:
            raise RuntimeError(
                f"Required component {required} was not retained"
            )


def validate_metrics(metrics_dir):
    for model in ["boltz1", "boltz2"]:
        for pose in range(5):
            path = metrics_dir / model / f"model_{pose}.json"

            if not path.is_file():
                raise RuntimeError(
                    f"Missing metrics file: {path}"
                )

            data = json.loads(path.read_text())
            rmsd = data.get("rmsd", {}).get(
                "assigned_scores", []
            )
            lddt = data.get("lddt_pli", {}).get(
                "assigned_scores", []
            )

            if not rmsd or not lddt:
                raise RuntimeError(
                    f"No valid ligand assignment in {path}. "
                    "Do not classify this pose as a failure."
                )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_id")
    parser.add_argument(
        "--root",
        type=Path,
        default=Path.home() / "boltz_visualization",
    )
    args = parser.parse_args()

    pdb_id = args.pdb_id.upper()
    root = args.root
    complex_dir = root / pdb_id
    config_path = complex_dir / "config.json"

    if not config_path.is_file():
        raise FileNotFoundError(config_path)

    config = json.loads(config_path.read_text())

    source_reference = (
        complex_dir / f"experimental_{pdb_id}.cif"
    )
    clean_reference_path = (
        complex_dir / f"experimental_{pdb_id}_clean.cif"
    )

    protein_chains = config["reference"]["protein_chains"]

    retained_nonpolymers = {
        ligand["ccd"]
        for ligand in config.get("ligands", [])
    }

    for cofactor in config.get("cofactors", []):
        if isinstance(cofactor, str):
            retained_nonpolymers.add(cofactor)
        else:
            retained_nonpolymers.add(cofactor["ccd"])

    for metal in config.get("metals", []):
        if isinstance(metal, str):
            retained_nonpolymers.add(metal)
        else:
            retained_nonpolymers.add(metal["ccd"])

    clean_reference(
        source_reference,
        clean_reference_path,
        protein_chains,
        retained_nonpolymers,
    )

    sdf_generator = root / "make_boltz_ligand_sdf.py"

    for model in ["boltz1", "boltz2"]:
        run([
            sys.executable,
            sdf_generator,
            "--pdb-id",
            pdb_id,
            "--model",
            model,
        ])

    metrics_dir = complex_dir / "metrics"

    for model in ["boltz1", "boltz2"]:
        (metrics_dir / model).mkdir(
            parents=True,
            exist_ok=True,
        )

        for pose in range(5):
            model_cif = (
                complex_dir
                / "production"
                / model
                / f"{pdb_id}_model_{pose}.cif"
            )
            ligand_sdf = (
                complex_dir
                / "production"
                / model
                / "ligand_sdf"
                / f"{pdb_id}_model_{pose}_ligand.sdf"
            )
            metric_output = (
                metrics_dir
                / model
                / f"model_{pose}.json"
            )

            comparison_command = [
                "ost",
                "compare-ligand-structures",
                "-m",
                model_cif,
                "-ml",
                ligand_sdf,
                "-r",
                clean_reference_path,
            ]

            reference_ligand_sdf = (
                complex_dir
                / f"experimental_{pdb_id}_ligand.sdf"
            )

            if reference_ligand_sdf.is_file():
                comparison_command.extend([
                    "-rl",
                    reference_ligand_sdf,
                ])

            comparison_command.extend([
                "-ft",
                "--lddt-pli",
                "--rmsd",
                "-fbs",
                "--full-results",
                "-o",
                metric_output,
            ])

            run(comparison_command)

    validate_metrics(metrics_dir)

    run([
        sys.executable,
        root / "summarize_boltz_results.py",
        "--pdb-id",
        pdb_id,
    ])

    print()
    print("=" * 60)
    print("Evaluation completed:", pdb_id)
    print(
        "Pose results:",
        complex_dir / f"{pdb_id}_pose_results.csv",
    )
    print(
        "Model summary:",
        complex_dir / f"{pdb_id}_model_summary.csv",
    )
    print("=" * 60)


if __name__ == "__main__":
    main()
