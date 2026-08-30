#!/usr/bin/env python3

import argparse
from collections import Counter
from pathlib import Path

import gemmi


WATER_NAMES = {"HOH", "DOD", "WAT"}


def values(block, tag):
    return [str(value).strip("'\"") for value in block.find_values(tag)]


def inspect(pdb_id, root):
    pdb_id = pdb_id.upper()
    directory = root / pdb_id
    cif_path = directory / f"experimental_{pdb_id}.cif"
    fasta_path = directory / f"{pdb_id}_all.fasta"

    print("\n" + "=" * 72)
    print("PDB:", pdb_id)
    print("=" * 72)

    if not cif_path.is_file() or cif_path.stat().st_size == 0:
        print("ERROR: missing or empty mmCIF:", cif_path)
        return

    if fasta_path.is_file():
        records = []
        header = None
        sequence = []

        for line in fasta_path.read_text().splitlines():
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence)))
                header = line[1:]
                sequence = []
            else:
                sequence.append(line.strip())

        if header is not None:
            records.append((header, "".join(sequence)))

        print("FASTA records:", len(records))
        for index, (header, sequence) in enumerate(records, start=1):
            print(f"  FASTA {index}: length={len(sequence)}")
            print(f"    {header}")
    else:
        print("WARNING: FASTA file not found")

    document = gemmi.cif.read_file(str(cif_path))
    block = document.sole_block()
    structure = gemmi.make_structure_from_block(block)

    chem_ids = values(block, "_chem_comp.id")
    chem_names = values(block, "_chem_comp.name")
    chem_formulas = values(block, "_chem_comp.formula")

    chemistry = {}
    for index, component_id in enumerate(chem_ids):
        chemistry[component_id] = {
            "name": chem_names[index] if index < len(chem_names) else "?",
            "formula": chem_formulas[index] if index < len(chem_formulas) else "?",
        }

    model = structure[0]
    observed = Counter()
    locations = {}

    print("\nCrystal chains:")

    for chain in model:
        protein_count = 0
        water_count = 0
        nonpolymers = []

        for residue in chain:
            name = residue.name.strip()

            if name in WATER_NAMES:
                water_count += 1
                continue

            if residue.het_flag == "A":
                protein_count += 1
            else:
                observed[name] += 1
                locations.setdefault(name, []).append(
                    f"{chain.name}:{residue.seqid}"
                )
                nonpolymers.append(f"{name}:{residue.seqid}")

        print(
            f"  chain={chain.name!r} "
            f"protein_residues={protein_count} "
            f"waters={water_count}"
        )

        if nonpolymers:
            print("    nonpolymers:", ", ".join(nonpolymers))

    print("\nObserved nonpolymer components:")

    if not observed:
        print("  None")
        return

    for component_id, count in sorted(observed.items()):
        information = chemistry.get(component_id, {})
        name = information.get("name", "?")
        formula = information.get("formula", "?")
        component_locations = ", ".join(locations.get(component_id, []))

        print(f"  CCD: {component_id}")
        print(f"    count: {count}")
        print(f"    locations: {component_locations}")
        print(f"    name: {name}")
        print(f"    formula: {formula}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_ids", nargs="+")
    args = parser.parse_args()

    root = Path.home() / "boltz_benchmark/reference"

    for pdb_id in args.pdb_ids:
        inspect(pdb_id, root)


if __name__ == "__main__":
    main()
