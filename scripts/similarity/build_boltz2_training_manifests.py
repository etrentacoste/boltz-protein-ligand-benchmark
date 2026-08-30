#!/usr/bin/env python3
"""Build versioned Boltz-2 training-system and PDB manifests."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--similarity-parquet", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--cutoff", default="2023-06-01")
    parser.add_argument(
        "--existing-pdb-manifest",
        type=Path,
        help="PDB IDs already present in the current structural reference.",
    )
    parser.add_argument(
        "--inclusive",
        action="store_true",
        help="Include entries released exactly on the cutoff date.",
    )
    return parser.parse_args()


def read_ids(path: Path) -> set[str]:
    return {
        line.strip().lower()
        for line in path.read_text().splitlines()
        if line.strip()
    }


def write_ids(path: Path, values: set[str] | list[str]) -> None:
    path.write_text("".join(f"{value}\n" for value in sorted(values)))


def main() -> None:
    args = parse_args()
    cutoff = pd.Timestamp(args.cutoff)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    frame = pd.read_parquet(
        args.similarity_parquet,
        columns=["target_system", "target_release_date"],
    )
    frame["target_system"] = frame["target_system"].astype("string").str.strip()
    frame["target_release_date"] = pd.to_datetime(
        frame["target_release_date"], errors="coerce"
    )
    frame = frame.dropna(subset=["target_system", "target_release_date"])
    frame = frame[frame["target_system"].ne("")]

    if args.inclusive:
        selected = frame[frame["target_release_date"].le(cutoff)].copy()
        cutoff_operator = "<="
    else:
        selected = frame[frame["target_release_date"].lt(cutoff)].copy()
        cutoff_operator = "<"

    # A system can occur in many query rows. Preserve its single earliest
    # recorded release date and write exactly one row per training system.
    systems = (
        selected.groupby("target_system", as_index=False)["target_release_date"]
        .min()
        .rename(columns={"target_system": "system_id"})
        .sort_values("system_id")
        .reset_index(drop=True)
    )
    systems["pdb_id"] = systems["system_id"].str.split("__", n=1).str[0].str.lower()
    invalid = ~systems["pdb_id"].str.fullmatch(r"[0-9a-z]{4}", na=False)
    if invalid.any():
        examples = systems.loc[invalid, "system_id"].head(10).tolist()
        raise RuntimeError(f"Invalid PLINDER system IDs: {examples}")

    system_ids = set(systems["system_id"])
    pdb_ids = set(systems["pdb_id"])
    systems.to_parquet(args.output_dir / "boltz2_training_systems.parquet", index=False)
    write_ids(args.output_dir / "boltz2_training_systems.txt", system_ids)
    write_ids(args.output_dir / "boltz2_training_pdb_ids.txt", pdb_ids)

    existing: set[str] = set()
    if args.existing_pdb_manifest is not None:
        existing = read_ids(args.existing_pdb_manifest)
    additional = pdb_ids - existing
    reused = pdb_ids & existing
    obsolete = existing - pdb_ids
    write_ids(args.output_dir / "boltz2_additional_pdb_ids.txt", additional)
    write_ids(args.output_dir / "boltz2_reused_pdb_ids.txt", reused)
    write_ids(args.output_dir / "existing_ids_not_in_boltz2.txt", obsolete)

    summary = pd.DataFrame(
        [
            {"item": "cutoff", "value": args.cutoff},
            {"item": "cutoff_operator", "value": cutoff_operator},
            {"item": "training_systems", "value": len(system_ids)},
            {"item": "training_pdb_ids", "value": len(pdb_ids)},
            {"item": "existing_pdb_ids", "value": len(existing)},
            {"item": "reused_pdb_ids", "value": len(reused)},
            {"item": "additional_pdb_ids", "value": len(additional)},
            {"item": "existing_ids_not_in_boltz2", "value": len(obsolete)},
            {
                "item": "minimum_release_date",
                "value": systems["target_release_date"].min(),
            },
            {
                "item": "maximum_release_date",
                "value": systems["target_release_date"].max(),
            },
        ]
    )
    summary.to_csv(args.output_dir / "manifest_summary.csv", index=False)

    print(f"Cutoff: target_release_date {cutoff_operator} {args.cutoff}")
    print(f"Training systems: {len(system_ids):,}")
    print(f"Training PDB IDs: {len(pdb_ids):,}")
    print(f"Existing PDB IDs: {len(existing):,}")
    print(f"Reused PDB IDs: {len(reused):,}")
    print(f"Additional PDB IDs: {len(additional):,}")
    print(f"Existing IDs outside Boltz-2 universe: {len(obsolete):,}")
    print(f"Output directory: {args.output_dir}")


if __name__ == "__main__":
    main()
