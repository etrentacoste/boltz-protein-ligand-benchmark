#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path

import gemmi


# Components that should not be selected as the principal ligand.
EXCLUDED_COMPONENTS = {
    # Water
    "HOH", "WAT", "DOD",

    # Common crystallisation solvents and buffers
    "EDO", "GOL", "PEG", "PGE", "PG4", "1PE",
    "MPD", "DMS", "DTT", "BME", "IPA", "EOH",
    "ACT", "ACY", "FMT", "MES", "TRS", "HEP",
    "BIS", "CAC", "CIT", "TAR", "MLA",

    # Common anions
    "SO4", "PO4", "NO3", "CO3", "SCN",
    "IOD", "BR", "CL", "F",

    # Common isolated metals and cations
    "NA", "K", "CA", "MG", "MN", "ZN", "FE",
    "CU", "CO", "NI", "CD", "HG", "AU",

    # Common glycans
    "NAG", "NDG", "MAN", "BMA", "FUC",
    "GAL", "GLC", "SIA",
}

# These may be functionally important, but normally should not be
# selected as the principal drug-like ligand.
COFACTOR_COMPONENTS = {
    "ATP", "ADP", "AMP",
    "GTP", "GDP", "GMP",
    "CTP", "CDP", "CMP",
    "UTP", "UDP", "UMP",
    "FAD", "FMN",
    "NAD", "NAP", "NAI", "NDP",
    "SAM", "SAH", "MTA",
    "PLP", "PMP",
    "HEM", "HEC",
    "COA", "ACO",
}

METAL_COMPONENTS = {
    "NA", "K", "CA", "MG", "MN", "ZN", "FE",
    "CU", "CO", "NI", "CD", "HG", "AU",
}


def clean_value(value):
    value = str(value).strip()

    if value in {"", ".", "?"}:
        return ""

    if (
        len(value) >= 2
        and value[0] == value[-1]
        and value[0] in {"'", '"'}
    ):
        value = value[1:-1]

    if value.startswith(";") and value.endswith(";"):
        value = value[1:-1].strip()

    return " ".join(value.split())


def values(block, tag):
    return [
        clean_value(value)
        for value in block.find_values(tag)
        if clean_value(value)
    ]


def first_value(block, tag):
    found = values(block, tag)
    return found[0] if found else ""


def parse_float(value):
    if not value:
        return None

    # Remove standard uncertainty, e.g. 1.234(5).
    value = re.sub(r"\([^)]*\)$", "", value)

    try:
        return float(value)
    except ValueError:
        return None


def formula_elements(formula):
    counts = {}

    for element, number in re.findall(
        r"([A-Z][a-z]?)(\d*)",
        formula,
    ):
        count = int(number) if number else 1
        counts[element] = counts.get(element, 0) + count

    return counts


def heavy_atom_count(formula):
    counts = formula_elements(formula)

    return sum(
        count
        for element, count in counts.items()
        if element not in {"H", "D"}
    )


def has_carbon(formula):
    return formula_elements(formula).get("C", 0) > 0


def split_chains(value):
    if not value:
        return []

    return [
        chain.strip()
        for chain in value.split(",")
        if chain.strip()
    ]


def chem_comp_table(block):
    # Read complete loop rows so that missing values represented
    # by "." or "?" do not shift the other columns.
    table = block.find(
        "_chem_comp.",
        [
            "id",
            "type",
            "name",
            "formula",
            "formula_weight",
        ],
    )

    result = {}

    for row in table:
        comp_id = clean_value(row[0]).upper()
        comp_type = clean_value(row[1])
        name = clean_value(row[2])
        formula = clean_value(row[3])
        weight = parse_float(clean_value(row[4]))

        if not comp_id:
            continue

        result[comp_id] = {
            "type": comp_type,
            "name": name,
            "formula": formula,
            "weight": weight,
            "heavy_atoms": (
                heavy_atom_count(formula)
                if formula
                else None
            ),
            "has_carbon": (
                has_carbon(formula)
                if formula
                else None
            ),
            "missing_formula": not bool(formula),
        }

    return result

def covalent_components(block):
    conn_types = values(
        block,
        "_struct_conn.conn_type_id",
    )
    partner_1 = values(
        block,
        "_struct_conn.ptnr1_label_comp_id",
    )
    partner_2 = values(
        block,
        "_struct_conn.ptnr2_label_comp_id",
    )

    if not conn_types:
        return set()

    if not (
        len(conn_types) == len(partner_1) == len(partner_2)
    ):
        return set()

    connected = set()

    for conn_type, comp_1, comp_2 in zip(
        conn_types,
        partner_1,
        partner_2,
    ):
        if conn_type.lower().startswith("covale"):
            connected.add(comp_1.upper())
            connected.add(comp_2.upper())

    return connected


def classify(path):
    block = gemmi.cif.read_file(str(path)).sole_block()

    pdb_id = first_value(block, "_entry.id").upper()
    title = first_value(block, "_struct.title")
    resolution = parse_float(
        first_value(block, "_refine.ls_d_res_high")
    )
    r_free = parse_float(
        first_value(block, "_refine.ls_R_factor_R_free")
    )

    revision_dates = values(
        block,
        "_pdbx_audit_revision_history.revision_date",
    )
    release_date = revision_dates[0] if revision_dates else ""

    entity_ids = values(block, "_entity_poly.entity_id")
    polymer_types = values(block, "_entity_poly.type")
    strand_ids = values(
        block,
        "_entity_poly.pdbx_strand_id",
    )
    mutations = values(
        block,
        "_entity_poly.pdbx_mutation",
    )

    protein_chains = []

    for entity_id, polymer_type, strands in zip(
        entity_ids,
        polymer_types,
        strand_ids,
    ):
        if "polypeptide" in polymer_type.lower():
            protein_chains.extend(split_chains(strands))

    protein_chains = sorted(set(protein_chains))

    mutation_values = [
        mutation
        for mutation in mutations
        if mutation.lower() not in {
            "no",
            "none",
            "not applicable",
        }
    ]

    nonpoly_ids = [
        value.upper()
        for value in values(
            block,
            "_pdbx_entity_nonpoly.comp_id",
        )
    ]

    components = chem_comp_table(block)
    covalent = covalent_components(block)

    candidate_ligands = []
    cofactors = []
    metals = []
    excluded_nonpolymers = []

    for comp_id in nonpoly_ids:
        component = components.get(comp_id)

        if component is None:
            excluded_nonpolymers.append(
                f"{comp_id}:missing_definition"
            )
            continue

        if comp_id in METAL_COMPONENTS:
            metals.append(comp_id)
            continue

        if comp_id in COFACTOR_COMPONENTS:
            cofactors.append(comp_id)
            continue

        if comp_id in EXCLUDED_COMPONENTS:
            excluded_nonpolymers.append(comp_id)
            continue

        if component["missing_formula"]:
            excluded_nonpolymers.append(
                f"{comp_id}:missing_formula"
            )
            continue

        if not component["has_carbon"]:
            excluded_nonpolymers.append(
                f"{comp_id}:no_carbon"
            )
            continue

        if component["heavy_atoms"] < 10:
            excluded_nonpolymers.append(
                f"{comp_id}:heavy_atoms_"
                f"{component['heavy_atoms']}"
            )
            continue

        weight = component["weight"]

        if weight is not None and not (120 <= weight <= 1200):
            excluded_nonpolymers.append(
                f"{comp_id}:mw_{weight:.1f}"
            )
            continue

        candidate_ligands.append(comp_id)

    # Preserve order while removing duplicates.
    candidate_ligands = list(
        dict.fromkeys(candidate_ligands)
    )
    cofactors = list(dict.fromkeys(cofactors))
    metals = list(dict.fromkeys(metals))

    reasons = []

    if not candidate_ligands:
        decision = "EXCLUDE"
        reasons.append("no_druglike_organic_ligand")

    elif len(candidate_ligands) > 1:
        decision = "MANUAL_REVIEW"
        reasons.append("multiple_candidate_ligands")

    else:
        target = candidate_ligands[0]

        if target in covalent:
            decision = "EXCLUDE"
            reasons.append("candidate_may_be_covalent")
        else:
            decision = "INCLUDE"

    if r_free is None:
        if decision == "INCLUDE":
            decision = "MANUAL_REVIEW"
        reasons.append("missing_r_free")

    elif r_free > 0.35:
        decision = "EXCLUDE"
        reasons.append("r_free_above_0.35")

    elif r_free > 0.28 and decision == "INCLUDE":
        decision = "RESERVE"
        reasons.append("r_free_above_0.28")

    if len(protein_chains) > 1:
        if decision == "INCLUDE":
            decision = "RESERVE"
        reasons.append(
            f"multiple_protein_chains_{len(protein_chains)}"
        )

    if mutation_values:
        if decision == "INCLUDE":
            decision = "RESERVE"
        reasons.append("reported_mutation")

    if cofactors:
        if decision == "INCLUDE":
            decision = "RESERVE"
        reasons.append("cofactor_present")

    if metals:
        if decision == "INCLUDE":
            decision = "RESERVE"
        reasons.append("metal_present")

    target = (
        candidate_ligands[0]
        if len(candidate_ligands) == 1
        else ""
    )

    target_data = components.get(target, {})

    return {
        "pdb_id": pdb_id,
        "release_date": release_date,
        "title": title,
        "resolution": (
            "" if resolution is None else resolution
        ),
        "r_free": "" if r_free is None else r_free,
        "protein_entity_count": sum(
            "polypeptide" in value.lower()
            for value in polymer_types
        ),
        "protein_chain_count": len(protein_chains),
        "protein_chains": ";".join(protein_chains),
        "reported_mutations": ";".join(mutation_values),
        "nonpolymer_components": ";".join(nonpoly_ids),
        "candidate_count": len(candidate_ligands),
        "candidate_ligands": ";".join(candidate_ligands),
        "target_ccd": target,
        "target_name": target_data.get("name", ""),
        "target_formula": target_data.get("formula", ""),
        "target_molecular_weight": (
            target_data.get("weight", "")
        ),
        "target_heavy_atoms": (
            target_data.get("heavy_atoms", "")
        ),
        "covalent_components": ";".join(
            sorted(covalent)
        ),
        "cofactors": ";".join(cofactors),
        "metals": ";".join(metals),
        "excluded_nonpolymers": ";".join(
            excluded_nonpolymers
        ),
        "decision": decision,
        "decision_reasons": ";".join(reasons),
    }


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--headers",
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
        input_rows = list(csv.DictReader(handle))

    rows = []

    for index, input_row in enumerate(
        input_rows,
        start=1,
    ):
        pdb_id = input_row["pdb_id"].upper()
        path = args.headers / f"{pdb_id}.cif"

        if not path.is_file():
            raise FileNotFoundError(path)

        row = classify(path)
        row["initial_rank"] = input_row.get("rank", index)
        rows.append(row)

        print(
            f"[{index:3d}/{len(input_rows)}] "
            f"{pdb_id}: {row['decision']} "
            f"{row['candidate_ligands']} "
            f"{row['decision_reasons']}"
        )

    fieldnames = [
        "initial_rank",
        "pdb_id",
        "release_date",
        "title",
        "resolution",
        "r_free",
        "protein_entity_count",
        "protein_chain_count",
        "protein_chains",
        "reported_mutations",
        "nonpolymer_components",
        "candidate_count",
        "candidate_ligands",
        "target_ccd",
        "target_name",
        "target_formula",
        "target_molecular_weight",
        "target_heavy_atoms",
        "covalent_components",
        "cofactors",
        "metals",
        "excluded_nonpolymers",
        "decision",
        "decision_reasons",
    ]

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
        writer.writerows(rows)

    print()
    print("Written:", args.output)
    print("Rows:", len(rows))

    for decision in [
        "INCLUDE",
        "RESERVE",
        "MANUAL_REVIEW",
        "EXCLUDE",
    ]:
        count = sum(
            row["decision"] == decision
            for row in rows
        )
        print(f"{decision}: {count}")


if __name__ == "__main__":
    main()
