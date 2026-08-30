#!/usr/bin/env python3

import argparse
import csv
import hashlib
from collections import Counter
from pathlib import Path

import gemmi


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
    result = []

    for value in block.find_values(tag):
        cleaned = clean_value(value)

        if cleaned:
            result.append(cleaned)

    return result


def sequence_signature(block):
    polymer_types = values(
        block,
        "_entity_poly.type",
    )
    sequences = values(
        block,
        "_entity_poly.pdbx_seq_one_letter_code_can",
    )

    protein_sequences = []

    if len(polymer_types) == len(sequences):
        for polymer_type, sequence in zip(
            polymer_types,
            sequences,
        ):
            if "polypeptide" not in polymer_type.lower():
                continue

            sequence = "".join(sequence.split()).upper()

            if sequence:
                protein_sequences.append(sequence)

    if not protein_sequences:
        return ""

    # Collapse identical copies in homomers.
    unique_sequences = sorted(set(protein_sequences))

    hashes = [
        hashlib.sha256(sequence.encode()).hexdigest()[:16]
        for sequence in unique_sequences
    ]

    return "SEQ:" + "+".join(hashes)


def target_signature(header_path):
    block = gemmi.cif.read_file(
        str(header_path)
    ).sole_block()

    db_names = values(block, "_struct_ref.db_name")
    accessions = values(
        block,
        "_struct_ref.pdbx_db_accession",
    )

    uniprot = []

    if len(db_names) == len(accessions):
        for db_name, accession in zip(
            db_names,
            accessions,
        ):
            normalized = db_name.upper()

            if (
                "UNP" in normalized
                or "UNIPROT" in normalized
            ):
                uniprot.append(accession.upper())

    uniprot = sorted(set(uniprot))

    if uniprot:
        return "UNIPROT:" + "+".join(uniprot)

    signature = sequence_signature(block)

    if signature:
        return signature

    return "PDB:" + header_path.stem.upper()


def number(value, default=999):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def reserve_score(row):
    score = 0
    chain_count = int(row["protein_chain_count"] or 0)
    r_free = number(row["r_free"])

    if row["metals"]:
        score += 4

    if row["cofactors"]:
        score += 4

    if chain_count == 2:
        score += 1
    elif chain_count > 2:
        score += 6

    if row["reported_mutations"]:
        score += 2

    if r_free > 0.30:
        score += 3
    elif r_free > 0.28:
        score += 1

    return score


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--classification",
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
    parser.add_argument(
        "--rejected-output",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--size",
        type=int,
        default=120,
    )
    parser.add_argument(
        "--max-per-target",
        type=int,
        default=2,
    )

    args = parser.parse_args()

    with args.classification.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    eligible = [
        row
        for row in rows
        if row["decision"] in {"INCLUDE", "RESERVE"}
    ]

    for row in eligible:
        pdb_id = row["pdb_id"]
        header = args.headers / f"{pdb_id}.cif"

        row["target_signature"] = target_signature(header)
        row["selection_priority"] = (
            0 if row["decision"] == "INCLUDE" else 1
        )
        row["reserve_score"] = (
            0
            if row["decision"] == "INCLUDE"
            else reserve_score(row)
        )

    eligible.sort(
        key=lambda row: (
            int(row["selection_priority"]),
            int(row["reserve_score"]),
            int(row["initial_rank"]),
        )
    )

    selected = []
    rejected = []
    signature_counts = Counter()

    for row in eligible:
        signature = row["target_signature"]

        if (
            signature_counts[signature]
            >= args.max_per_target
        ):
            copied = dict(row)
            copied["pool_status"] = (
                "REJECTED_TARGET_REDUNDANCY"
            )
            copied["pool_rejection_reason"] = (
                f"more_than_{args.max_per_target}_"
                "entries_for_target"
            )
            rejected.append(copied)
            continue

        if len(selected) >= args.size:
            copied = dict(row)
            copied["pool_status"] = "NOT_SELECTED"
            copied["pool_rejection_reason"] = (
                "pool_size_reached"
            )
            rejected.append(copied)
            continue

        signature_counts[signature] += 1

        copied = dict(row)
        copied["pool_status"] = "SELECTED"
        copied["pool_rejection_reason"] = ""
        copied["target_instance_in_pool"] = (
            signature_counts[signature]
        )
        selected.append(copied)

    if len(selected) < args.size:
        raise RuntimeError(
            f"Only {len(selected)} non-redundant entries "
            f"could be selected; requested {args.size}"
        )

    all_fields = []

    for row in selected + rejected:
        for field in row:
            if field not in all_fields:
                all_fields.append(field)

    args.output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=all_fields,
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(selected)

    with args.rejected_output.open(
        "w",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=all_fields,
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rejected)

    print("Eligible INCLUDE/RESERVE:", len(eligible))
    print("Selected:", len(selected))
    print("Rejected/not selected:", len(rejected))
    print(
        "Unique target signatures:",
        len(signature_counts),
    )

    repeated = {
        signature: count
        for signature, count in signature_counts.items()
        if count > 1
    }

    print(
        "Targets represented twice:",
        len(repeated),
    )
    print("Written:", args.output)
    print("Rejected log:", args.rejected_output)


if __name__ == "__main__":
    main()
