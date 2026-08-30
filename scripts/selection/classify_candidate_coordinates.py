#!/usr/bin/env python3

import argparse
import csv
import json
from pathlib import Path

import gemmi
import numpy as np


WATER_NAMES = {"HOH", "WAT", "DOD"}

STANDARD_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
    "SEC", "PYL",
}


def split_values(value):
    return {
        item.strip().upper()
        for item in str(value).split(";")
        if item.strip()
    }


def atom_is_hydrogen(atom):
    return atom.element.name.upper() in {"H", "D"}


def atom_altloc(atom):
    value = str(atom.altloc).strip()

    if value in {"", ".", "?", "\x00"}:
        return ""

    return value


def residue_is_protein(residue):
    if residue.name.upper() in STANDARD_AMINO_ACIDS:
        return True

    try:
        return (
            residue.entity_type
            == gemmi.EntityType.Polymer
            and residue.name.upper() not in WATER_NAMES
        )
    except Exception:
        return (
            residue.het_flag != "H"
            and residue.name.upper() not in WATER_NAMES
        )


def atom_coordinates(atoms):
    return np.asarray(
        [
            [atom.pos.x, atom.pos.y, atom.pos.z]
            for atom in atoms
        ],
        dtype=float,
    )


def minimum_distance(coords_a, coords_b):
    if len(coords_a) == 0 or len(coords_b) == 0:
        return None

    minimum_squared = float("inf")

    # Work one atom at a time to avoid a very large
    # ligand x protein x 3 matrix.
    for coordinate in coords_a:
        differences = coords_b - coordinate
        squared = np.einsum(
            "ij,ij->i",
            differences,
            differences,
        )

        current = float(np.min(squared))

        if current < minimum_squared:
            minimum_squared = current

    return minimum_squared ** 0.5


def ligand_protein_contacts(
    ligand_coords,
    protein_coords,
    protein_metadata,
    cutoff=5.0,
):
    if (
        len(ligand_coords) == 0
        or len(protein_coords) == 0
    ):
        return None, set(), set()

    cutoff_squared = cutoff * cutoff
    minimum_squared = float("inf")
    contacting_residues = set()
    contacting_chains = set()

    for coordinate in ligand_coords:
        differences = protein_coords - coordinate
        squared = np.einsum(
            "ij,ij->i",
            differences,
            differences,
        )

        current = float(np.min(squared))

        if current < minimum_squared:
            minimum_squared = current

        indices = np.where(squared <= cutoff_squared)[0]

        for index in indices:
            chain_name, residue_name, residue_number = (
                protein_metadata[index]
            )
            contacting_chains.add(chain_name)
            contacting_residues.add(
                (
                    chain_name,
                    residue_name,
                    residue_number,
                )
            )

    return (
        minimum_squared ** 0.5,
        contacting_residues,
        contacting_chains,
    )


def collect_component_atoms(model, component_ids):
    result = {}

    for component_id in component_ids:
        result[component_id] = []

    for chain in model:
        for residue in chain:
            name = residue.name.upper()

            if name not in component_ids:
                continue

            atoms = [
                atom
                for atom in residue
                if not atom_is_hydrogen(atom)
            ]

            if atoms:
                result[name].extend(atoms)

    return result


def analyze_structure(path, target_ccd, metals, cofactors):
    structure = gemmi.read_structure(str(path))

    if len(structure) == 0:
        raise RuntimeError("Structure has no models")

    model = structure[0]

    protein_atoms = []
    protein_metadata = []
    ligand_instances = []

    for chain in model:
        for residue in chain:
            heavy_atoms = [
                atom
                for atom in residue
                if not atom_is_hydrogen(atom)
            ]

            if not heavy_atoms:
                continue

            if residue_is_protein(residue):
                residue_number = str(residue.seqid)

                for atom in heavy_atoms:
                    protein_atoms.append(atom)
                    protein_metadata.append(
                        (
                            chain.name,
                            residue.name,
                            residue_number,
                        )
                    )

            if residue.name.upper() == target_ccd:
                ligand_instances.append(
                    (
                        chain.name,
                        residue,
                        heavy_atoms,
                    )
                )

    protein_coords = atom_coordinates(protein_atoms)

    nearby_component_ids = metals | cofactors
    component_atoms = collect_component_atoms(
        model,
        nearby_component_ids,
    )

    component_coords = {
        component_id: atom_coordinates(atoms)
        for component_id, atoms in component_atoms.items()
        if atoms
    }

    instance_results = []

    for chain_name, residue, atoms in ligand_instances:
        ligand_coords = atom_coordinates(atoms)

        (
            protein_minimum,
            contact_residues,
            contact_chains,
        ) = ligand_protein_contacts(
            ligand_coords,
            protein_coords,
            protein_metadata,
        )

        occupancies = [
            float(atom.occ)
            for atom in atoms
        ]
        b_factors = [
            float(atom.b_iso)
            for atom in atoms
        ]
        altlocs = sorted({
            atom_altloc(atom)
            for atom in atoms
            if atom_altloc(atom)
        })

        nearby_distances = {}

        for component_id, coordinates in component_coords.items():
            distance = minimum_distance(
                ligand_coords,
                coordinates,
            )

            if distance is not None:
                nearby_distances[component_id] = distance

        result = {
            "chain": chain_name,
            "residue_number": str(residue.seqid),
            "heavy_atoms": len(atoms),
            "mean_occupancy": float(np.mean(occupancies)),
            "minimum_occupancy": float(np.min(occupancies)),
            "mean_b_factor": float(np.mean(b_factors)),
            "altlocs": altlocs,
            "protein_minimum_distance": protein_minimum,
            "contact_residue_count": len(contact_residues),
            "contact_chain_count": len(contact_chains),
            "contact_chains": sorted(contact_chains),
            "contact_residues": [
                f"{chain}:{name}:{number}"
                for chain, name, number
                in sorted(contact_residues)
            ],
            "component_distances": nearby_distances,
        }

        # Prefer an instance with a well-defined and extensive
        # protein-binding environment.
        result["_ranking"] = (
            -result["contact_residue_count"],
            -result["mean_occupancy"],
            result["mean_b_factor"],
            chain_name,
            str(residue.seqid),
        )

        instance_results.append(result)

    instance_results.sort(
        key=lambda result: result["_ranking"]
    )

    for result in instance_results:
        del result["_ranking"]

    return instance_results


def classify(row, instances):
    reasons = []
    notes = []

    if not instances:
        return {
            "coordinate_decision": "EXCLUDE",
            "coordinate_reasons": "target_ligand_not_found",
        }

    best = instances[0]

    mean_occupancy = best["mean_occupancy"]
    minimum_occupancy = best["minimum_occupancy"]
    contact_count = best["contact_residue_count"]
    contact_chain_count = best["contact_chain_count"]
    protein_minimum = best["protein_minimum_distance"]

    decision = "INCLUDE"

    if protein_minimum is None or contact_count < 3:
        decision = "EXCLUDE"
        reasons.append("insufficient_protein_contacts")

    elif contact_count < 5:
        decision = "RESERVE"
        reasons.append("few_protein_contact_residues")

    if mean_occupancy < 0.70:
        decision = "EXCLUDE"
        reasons.append("mean_occupancy_below_0.70")

    elif mean_occupancy < 0.90:
        if decision == "INCLUDE":
            decision = "RESERVE"
        reasons.append("mean_occupancy_below_0.90")

    if minimum_occupancy < 0.50:
        decision = "EXCLUDE"
        reasons.append("minimum_occupancy_below_0.50")

    elif minimum_occupancy < 0.90:
        if decision == "INCLUDE":
            decision = "RESERVE"
        reasons.append("some_atoms_occupancy_below_0.90")

    if best["altlocs"]:
        if decision == "INCLUDE":
            decision = "RESERVE"
        reasons.append("ligand_alternate_conformation")

    if contact_chain_count > 1:
        if decision == "INCLUDE":
            decision = "RESERVE"
        reasons.append(
            f"binding_site_uses_{contact_chain_count}_chains"
        )

    if len(instances) > 1:
        if decision == "INCLUDE":
            decision = "RESERVE"
        reasons.append(
            f"multiple_ligand_instances_{len(instances)}"
        )

    component_distances = best["component_distances"]

    metals = split_values(row.get("metals", ""))
    cofactors = split_values(row.get("cofactors", ""))

    nearby_metals = {
        component: distance
        for component, distance in component_distances.items()
        if component in metals and distance <= 5.0
    }

    nearby_cofactors = {
        component: distance
        for component, distance in component_distances.items()
        if component in cofactors and distance <= 5.0
    }

    if nearby_metals:
        if decision == "INCLUDE":
            decision = "RESERVE"
        reasons.append("metal_within_5A")

    if nearby_cofactors:
        if decision == "INCLUDE":
            decision = "RESERVE"
        reasons.append("cofactor_within_5A")

    try:
        r_free = float(row["r_free"])
    except (TypeError, ValueError):
        r_free = None

    if r_free is None:
        if decision == "INCLUDE":
            decision = "MANUAL_REVIEW"
        reasons.append("missing_r_free")

    elif r_free > 0.35:
        decision = "EXCLUDE"
        reasons.append("r_free_above_0.35")

    elif r_free > 0.28:
        if decision == "INCLUDE":
            decision = "RESERVE"
        reasons.append("r_free_above_0.28")

    if row.get("reported_mutations"):
        if decision == "INCLUDE":
            decision = "RESERVE"
        reasons.append("reported_mutation")

    # Metal/cofactor present in the structure but distant from the
    # selected ligand should not by itself force RESERVE.
    if metals and not nearby_metals:
        notes.append("all_metals_farther_than_5A")

    if cofactors and not nearby_cofactors:
        notes.append("all_cofactors_farther_than_5A")

    if not reasons:
        reasons.append("clean_coordinate_environment")

    return {
        "coordinate_decision": decision,
        "coordinate_reasons": ";".join(reasons),
        "coordinate_notes": ";".join(notes),
        "selected_ligand_chain": best["chain"],
        "selected_ligand_residue": best["residue_number"],
        "ligand_instance_count": len(instances),
        "ligand_heavy_atoms_observed": best["heavy_atoms"],
        "ligand_mean_occupancy": (
            f"{best['mean_occupancy']:.4f}"
        ),
        "ligand_minimum_occupancy": (
            f"{best['minimum_occupancy']:.4f}"
        ),
        "ligand_mean_b_factor": (
            f"{best['mean_b_factor']:.3f}"
        ),
        "ligand_altlocs": ";".join(best["altlocs"]),
        "protein_minimum_distance": (
            ""
            if best["protein_minimum_distance"] is None
            else f"{best['protein_minimum_distance']:.3f}"
        ),
        "contact_residue_count": (
            best["contact_residue_count"]
        ),
        "contact_chain_count": (
            best["contact_chain_count"]
        ),
        "contact_chains": ";".join(
            best["contact_chains"]
        ),
        "contact_residues": ";".join(
            best["contact_residues"]
        ),
        "nearby_metals": json.dumps(
            nearby_metals,
            sort_keys=True,
        ),
        "nearby_cofactors": json.dumps(
            nearby_cofactors,
            sort_keys=True,
        ),
        "all_ligand_instances_json": json.dumps(
            instances,
            sort_keys=True,
        ),
    }


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--cif-dir",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
    )

    args = parser.parse_args()

    with args.input.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    output_rows = []

    for index, row in enumerate(rows, start=1):
        pdb_id = row["pdb_id"].upper()
        target_ccd = row["target_ccd"].upper()
        path = args.cif_dir / f"{pdb_id}.cif"

        metals = split_values(row.get("metals", ""))
        cofactors = split_values(
            row.get("cofactors", "")
        )

        instances = analyze_structure(
            path,
            target_ccd,
            metals,
            cofactors,
        )

        result = classify(row, instances)
        combined = dict(row)
        combined.update(result)
        output_rows.append(combined)

        print(
            f"[{index:3d}/{len(rows)}] "
            f"{pdb_id} {target_ccd}: "
            f"{result['coordinate_decision']} "
            f"{result['coordinate_reasons']}"
        )

    fieldnames = []

    for row in output_rows:
        for field in row:
            if field not in fieldnames:
                fieldnames.append(field)

    args.output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(output_rows)

    print()
    print("Written:", args.output)
    print("Rows:", len(output_rows))

    for decision in [
        "INCLUDE",
        "RESERVE",
        "MANUAL_REVIEW",
        "EXCLUDE",
    ]:
        count = sum(
            row["coordinate_decision"] == decision
            for row in output_rows
        )
        print(f"{decision}: {count}")


if __name__ == "__main__":
    main()
