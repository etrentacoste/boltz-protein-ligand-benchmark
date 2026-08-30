#!/usr/bin/env python3
"""Summarize per-query PLINDER target manifests and build their union."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def read_ids(path: Path) -> set[str]:
    return {line.strip() for line in path.read_text().splitlines() if line.strip()}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--query-manifest", required=True, type=Path)
    parser.add_argument("--enumeration-dir", required=True, type=Path)
    parser.add_argument("--systems-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    query_ids = sorted(read_ids(args.query_manifest))
    union: set[str] = set()
    rows: list[dict[str, object]] = []
    missing_manifests: list[str] = []

    for pdb_id in query_ids:
        manifest = (
            args.enumeration_dir / pdb_id / "independent_target_systems.txt"
        )
        if not manifest.is_file():
            missing_manifests.append(pdb_id)
            continue
        systems = read_ids(manifest)
        union.update(systems)
        present = sum((args.systems_dir / system_id).is_dir() for system_id in systems)
        rows.append(
            {
                "query_pdb_id": pdb_id,
                "candidate_systems": len(systems),
                "already_materialized": present,
                "missing_systems": len(systems) - present,
            }
        )

    present_union = {
        system_id for system_id in union if (args.systems_dir / system_id).is_dir()
    }
    missing_union = union - present_union

    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).sort_values("query_pdb_id").to_csv(
        args.output_dir / "per_query_target_summary.csv", index=False
    )
    (args.output_dir / "all_candidate_systems.txt").write_text(
        "".join(f"{system_id}\n" for system_id in sorted(union))
    )
    (args.output_dir / "missing_candidate_systems.txt").write_text(
        "".join(f"{system_id}\n" for system_id in sorted(missing_union))
    )

    print(f"Queries requested: {len(query_ids):,}")
    print(f"Queries mapped: {len(rows):,}")
    print(f"Missing query manifests: {len(missing_manifests):,}")
    if missing_manifests:
        print("Missing:", ", ".join(missing_manifests))
    print(f"Unique candidate systems: {len(union):,}")
    print(f"Already materialized: {len(present_union):,}")
    print(f"Additional systems required: {len(missing_union):,}")
    print(f"Output directory: {args.output_dir}")


if __name__ == "__main__":
    main()
