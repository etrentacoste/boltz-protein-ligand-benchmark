import argparse
import json
from pathlib import Path

import yaml


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=Path)
    parser.add_argument(
        "--benchmark-root",
        type=Path,
        default=Path.home() / "boltz_benchmark",
    )
    args = parser.parse_args()

    config = json.loads(args.config.read_text())
    pdb_id = config["pdb_id"].upper()

    sequences_with_msa = []
    sequences_without_msa = []

    for protein in config["proteins"]:
        with_msa = {
            "protein": {
                "id": protein["id"],
                "sequence": protein["sequence"],
                "msa": protein["msa"],
            }
        }

        without_msa = {
            "protein": {
                "id": protein["id"],
                "sequence": protein["sequence"],
            }
        }

        sequences_with_msa.append(with_msa)
        sequences_without_msa.append(without_msa)

    binder_ids = []

    for ligand in config["ligands"]:
        ligand_data = {
            "id": ligand["id"],
        }

        has_smiles = bool(ligand.get("smiles"))
        has_ccd = bool(ligand.get("ccd"))

        if not has_smiles and not has_ccd:
            raise ValueError(
                "Each ligand must provide at least one of "
                f"'smiles' or 'ccd': {ligand}"
            )

        # Prefer SMILES when both are retained in the benchmark
        # metadata. Use CCD for cofactors and ions without SMILES.
        if has_smiles:
            ligand_data["smiles"] = ligand["smiles"]
        else:
            ligand_data["ccd"] = ligand["ccd"]

        sequences_with_msa.append(
            {"ligand": dict(ligand_data)}
        )
        sequences_without_msa.append(
            {"ligand": dict(ligand_data)}
        )

        if ligand.get("affinity_binder", False):
            binder_ids.append(ligand["id"])

    if len(binder_ids) != 1:
        raise ValueError(
            "Exactly one affinity binder must be specified; "
            f"found {binder_ids}"
        )

    boltz1 = {
        "version": 1,
        "sequences": sequences_with_msa,
    }

    boltz2 = {
        "version": 1,
        "sequences": sequences_with_msa,
        "properties": [
            {
                "affinity": {
                    "binder": binder_ids[0],
                }
            }
        ],
    }

    msa_generation = {
        "version": 1,
        "sequences": sequences_without_msa,
    }

    outputs = {
        (
            args.benchmark_root
            / f"inputs/boltz1/{pdb_id}.yaml"
        ): boltz1,
        (
            args.benchmark_root
            / f"inputs/boltz2/{pdb_id}.yaml"
        ): boltz2,
        (
            args.benchmark_root
            / f"inputs/msa_generation/{pdb_id}.yaml"
        ): msa_generation,
    }

    for output, data in outputs.items():
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            yaml.safe_dump(
                data,
                sort_keys=False,
                width=10000,
            )
        )
        print("Written:", output)


if __name__ == "__main__":
    main()
