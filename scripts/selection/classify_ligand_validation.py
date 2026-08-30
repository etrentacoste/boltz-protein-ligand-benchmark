#!/usr/bin/env python3

import argparse
import csv
import gzip
import re
import xml.etree.ElementTree as ET
from pathlib import Path


def parse_float(value):
    if value in {None, "", ".", "?"}:
        return None

    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def parse_residue_identifier(value):
    value = str(value).strip()

    match = re.search(
        r"(-?\d+)\s*([A-Za-z]?)",
        value,
    )

    if not match:
        return value, ""

    return match.group(1), match.group(2)


def validation_records(report_path, target_ccd):
    with gzip.open(report_path, "rb") as handle:
        root = ET.parse(handle).getroot()

    records = []

    for element in root.iter():
        tag = element.tag.split("}")[-1]

        if tag != "ModelledSubgroup":
            continue

        attributes = dict(element.attrib)

        if (
            attributes.get("resname", "").upper()
            != target_ccd.upper()
        ):
            continue

        records.append(attributes)

    return records


def choose_record(
    records,
    selected_chain,
    selected_residue,
):
    residue_number, insertion_code = (
        parse_residue_identifier(selected_residue)
    )

    exact = [
        record
        for record in records
        if (
            record.get("chain", "") == selected_chain
            and record.get("resnum", "") == residue_number
            and record.get("icode", "") == insertion_code
        )
    ]

    if exact:
        candidates = exact
        match_quality = "exact_chain_residue"
    elif len(records) == 1:
        candidates = records
        match_quality = "single_target_record_fallback"
    else:
        same_chain = [
            record
            for record in records
            if record.get("chain", "") == selected_chain
        ]

        if same_chain:
            candidates = same_chain
            match_quality = "chain_only_fallback"
        else:
            candidates = records
            match_quality = "best_available_fallback"

    def ranking(record):
        empty_alt = (
            1
            if record.get("altcode", "") in {"", ".", "?"}
            else 0
        )
        model_one = (
            1
            if record.get("model", "1") == "1"
            else 0
        )
        rscc = parse_float(record.get("rscc"))
        occupancy = parse_float(record.get("avgoccu"))

        return (
            model_one,
            empty_alt,
            -1 if rscc is None else rscc,
            -1 if occupancy is None else occupancy,
        )

    return max(candidates, key=ranking), match_quality


def classify(row, report_path):
    target_ccd = row["target_ccd"].upper()
    coordinate_decision = row["coordinate_decision"]

    records = validation_records(
        report_path,
        target_ccd,
    )

    if not records:
        if coordinate_decision == "EXCLUDE":
            decision = "EXCLUDE"
        else:
            decision = "MANUAL_REVIEW"

        return {
            "validation_decision": decision,
            "validation_reasons": (
                "no_target_validation_record"
            ),
            "validation_record_count": 0,
        }

    record, match_quality = choose_record(
        records,
        row["selected_ligand_chain"],
        row["selected_ligand_residue"],
    )

    rscc = parse_float(record.get("rscc"))
    rsr = parse_float(record.get("rsr"))
    occupancy = parse_float(record.get("avgoccu"))
    ediam = parse_float(record.get("EDIAm"))
    opia = parse_float(record.get("OPIA"))
    owab = parse_float(record.get("owab"))
    natoms_eds = parse_float(record.get("NatomsEDS"))

    bond_rmsz = parse_float(
        record.get("mogul_bonds_rmsz")
    )
    angle_rmsz = parse_float(
        record.get("mogul_angles_rmsz")
    )

    reasons = []

    if coordinate_decision == "EXCLUDE":
        decision = "EXCLUDE"
        reasons.append("coordinate_screen_exclude")

    elif coordinate_decision == "RESERVE":
        decision = "RESERVE"
        reasons.append("coordinate_screen_reserve")

    else:
        decision = "INCLUDE"

    if rscc is None:
        if decision == "INCLUDE":
            decision = "MANUAL_REVIEW"
        reasons.append("missing_rscc")

    elif rscc < 0.70:
        decision = "EXCLUDE"
        reasons.append("rscc_below_0.70")

    elif rscc < 0.85:
        if decision == "INCLUDE":
            decision = "RESERVE"
        reasons.append("rscc_below_0.85")

    if rsr is None:
        if decision == "INCLUDE":
            decision = "MANUAL_REVIEW"
        reasons.append("missing_rsr")

    elif rsr > 0.40:
        decision = "EXCLUDE"
        reasons.append("rsr_above_0.40")

    elif rsr > 0.30:
        if decision == "INCLUDE":
            decision = "RESERVE"
        reasons.append("rsr_above_0.30")

    if occupancy is None:
        if decision == "INCLUDE":
            decision = "MANUAL_REVIEW"
        reasons.append("missing_validation_occupancy")

    elif occupancy < 0.70:
        decision = "EXCLUDE"
        reasons.append(
            "validation_occupancy_below_0.70"
        )

    elif occupancy < 0.90:
        if decision == "INCLUDE":
            decision = "RESERVE"
        reasons.append(
            "validation_occupancy_below_0.90"
        )

    if not reasons:
        reasons.append("validation_quality_pass")

    observed_heavy_atoms = parse_float(
        row.get("ligand_heavy_atoms_observed")
    )

    eds_atom_coverage = None

    if (
        natoms_eds is not None
        and observed_heavy_atoms not in {None, 0}
    ):
        eds_atom_coverage = (
            natoms_eds / observed_heavy_atoms
        )

    return {
        "validation_decision": decision,
        "validation_reasons": ";".join(reasons),
        "validation_record_count": len(records),
        "validation_match_quality": match_quality,
        "validation_chain": record.get("chain", ""),
        "validation_resnum": record.get("resnum", ""),
        "validation_icode": record.get("icode", ""),
        "validation_altcode": record.get("altcode", ""),
        "validation_model": record.get("model", ""),
        "ligand_rscc": "" if rscc is None else rscc,
        "ligand_rsr": "" if rsr is None else rsr,
        "validation_avg_occupancy": (
            "" if occupancy is None else occupancy
        ),
        "ligand_ediam": "" if ediam is None else ediam,
        "ligand_opia": "" if opia is None else opia,
        "ligand_owab": "" if owab is None else owab,
        "ligand_natoms_eds": (
            "" if natoms_eds is None else natoms_eds
        ),
        "ligand_eds_atom_coverage": (
            ""
            if eds_atom_coverage is None
            else f"{eds_atom_coverage:.4f}"
        ),
        "mogul_bonds_rmsz": (
            "" if bond_rmsz is None else bond_rmsz
        ),
        "mogul_angles_rmsz": (
            "" if angle_rmsz is None else angle_rmsz
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
        "--report-dir",
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
        report = (
            args.report_dir
            / f"{pdb_id}_validation.xml.gz"
        )

        if not report.is_file():
            result = {
                "validation_decision": "MANUAL_REVIEW",
                "validation_reasons": (
                    "validation_report_missing"
                ),
                "validation_record_count": 0,
            }
        else:
            result = classify(row, report)

        combined = dict(row)
        combined.update(result)
        output_rows.append(combined)

        print(
            f"[{index:3d}/{len(rows)}] "
            f"{pdb_id} {row['target_ccd']}: "
            f"{result['validation_decision']} "
            f"{result['validation_reasons']}"
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
            row["validation_decision"] == decision
            for row in output_rows
        )
        print(f"{decision}: {count}")


if __name__ == "__main__":
    main()
