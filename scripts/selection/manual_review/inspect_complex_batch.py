import argparse
import shlex
from collections import Counter
from pathlib import Path

import gemmi


WATER = {"HOH", "WAT", "DOD"}


def read_fasta(path):
    records = []
    header = None
    sequence = []

    for line in path.read_text().splitlines():
        line = line.strip()

        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(sequence)))
            header = line[1:]
            sequence = []
        elif line:
            sequence.append(line)

    if header is not None:
        records.append((header, "".join(sequence)))

    return records


def read_smiles(path):
    candidates = []

    for line in path.read_text().splitlines():
        if "SMILES" not in line:
            continue

        try:
            fields = shlex.split(line)
        except ValueError:
            continue

        if len(fields) < 5:
            continue

        if fields[0] != path.stem:
            continue

        descriptor_type = fields[1]
        program = fields[2]
        descriptor = fields[-1]

        candidates.append(
            (descriptor_type, program, descriptor)
        )

    preferred = [
        item for item in candidates
        if (
            item[0] == "SMILES_CANONICAL"
            and "OpenEye" in item[1]
        )
    ]

    if preferred:
        return preferred[0][2], candidates

    canonical = [
        item for item in candidates
        if item[0] == "SMILES_CANONICAL"
    ]

    if canonical:
        return canonical[0][2], candidates

    return None, candidates


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root",
        type=Path,
        default=Path.home() / "boltz_benchmark/reference",
    )
    args = parser.parse_args()

    systems = {
        "9TKU": "A1JWM",
        "9LOK": "A1L7F",
        "9RQV": "A1JIG",
        "9RMF": "A1JPF",
        "9Y56": "A1CSI",
    }

    for pdb_id, ccd in systems.items():
        directory = args.root / pdb_id
        fasta = directory / f"{pdb_id}_all.fasta"
        cif = directory / f"experimental_{pdb_id}.cif"
        ccd_file = directory / f"{ccd}.cif"

        records = read_fasta(fasta)
        smiles, candidates = read_smiles(ccd_file)
        structure = gemmi.read_structure(str(cif))
        model = structure[0]

        print()
        print("=" * 72)
        print("PDB:", pdb_id)
        print("Target CCD:", ccd)
        print("FASTA records:", len(records))

        for index, (header, sequence) in enumerate(records, 1):
            print(
                f"  FASTA {index}: length={len(sequence)} "
                f"header={header}"
            )

        print("Selected SMILES:", smiles)
        print("Crystal chains:")

        for chain in model:
            nonstandard = Counter(
                residue.name
                for residue in chain
                if (
                    residue.het_flag == "H"
                    and residue.name not in WATER
                )
            )

            target_instances = [
                str(residue.seqid)
                for residue in chain
                if residue.name == ccd
            ]

            protein_residues = sum(
                residue.het_flag != "H"
                for residue in chain
            )
            waters = sum(
                residue.name in WATER
                for residue in chain
            )

            print(
                f"  chain={chain.name} "
                f"protein_residues={protein_residues} "
                f"waters={waters} "
                f"target_instances={target_instances} "
                f"other_nonpolymer={dict(nonstandard)}"
            )


if __name__ == "__main__":
    main()
