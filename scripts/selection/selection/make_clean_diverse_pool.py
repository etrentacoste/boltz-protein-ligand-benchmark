#!/usr/bin/env python3

import argparse
import csv
import importlib.util
from collections import Counter
from pathlib import Path


def load_signature_function(script_path):
    spec = importlib.util.spec_from_file_location(
        "diverse_pool_helpers",
        script_path,
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.target_signature


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
        "--signature-script",
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
        "--max-per-target",
        type=int,
        default=2,
    )

    args = parser.parse_args()

    target_signature = load_signature_function(
        args.signature_script
    )

    with args.classification.open(newline="") as handle:
        rows = [
            row
            for row in csv.DictReader(handle)
            if row["decision"] == "INCLUDE"
        ]

    rows.sort(
        key=lambda row: int(row["initial_rank"])
    )

    signature_counts = Counter()
    selected = []
    rejected = []

    for row in rows:
        pdb_id = row["pdb_id"].upper()
        header = args.headers / f"{pdb_id}.cif"
        signature = target_signature(header)

        copied = dict(row)
        copied["target_signature"] = signature

        if (
            signature_counts[signature]
            >= args.max_per_target
        ):
            copied["diversity_decision"] = (
                "REJECT_REDUNDANCY"
            )
            copied["diversity_reason"] = (
                f"more_than_{args.max_per_target}_"
                "entries_for_target"
            )
            rejected.append(copied)
            continue

        signature_counts[signature] += 1

        copied["diversity_decision"] = "SELECT"
        copied["diversity_reason"] = ""
        copied["target_instance_number"] = (
            signature_counts[signature]
        )
        selected.append(copied)

    fieldnames = []

    for row in selected + rejected:
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
            fieldnames=fieldnames,
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rejected)

    print("Header INCLUDE:", len(rows))
    print(
        "Selected after diversity:",
        len(selected),
    )
    print(
        "Rejected for target redundancy:",
        len(rejected),
    )
    print(
        "Unique target signatures:",
        len(signature_counts),
    )
    print(
        "Maximum entries per target:",
        max(signature_counts.values()),
    )
    print("Written:", args.output)
    print("Rejected log:", args.rejected_output)


if __name__ == "__main__":
    main()
