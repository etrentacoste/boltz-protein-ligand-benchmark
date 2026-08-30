#!/usr/bin/env python3
"""Calculate the RDKit-2024 SuCOS-shape metric for one query PDB."""

from __future__ import annotations

import argparse
from pathlib import Path

import calculate_validate_sucos_shape as sucos
import numpy as np
import pandas as pd
from rdkit import Chem


def read_manifest(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--query-systems-dir", required=True, type=Path)
    parser.add_argument("--target-systems-dir", required=True, type=Path)
    parser.add_argument("--query-manifest", required=True, type=Path)
    parser.add_argument("--target-manifest", required=True, type=Path)
    parser.add_argument("--output-parquet", required=True, type=Path)
    return parser.parse_args()


def read_ligands(directory: Path) -> tuple[list[Chem.Mol], int]:
    molecules = []
    failures = 0
    for path in sorted((directory / "ligand_files").glob("*.sdf")):
        molecule = sucos.read_molecule(path)
        if molecule is None:
            failures += 1
        else:
            molecules.append(molecule)
    return molecules, failures


def main() -> None:
    args = parse_args()
    query_ids = read_manifest(args.query_manifest)
    target_ids = read_manifest(args.target_manifest)
    query_molecules = {}
    query_read_failures = 0
    for system_id in query_ids:
        molecules, failures = read_ligands(args.query_systems_dir / system_id)
        query_molecules[system_id] = molecules
        query_read_failures += failures
    missing = [system_id for system_id, mols in query_molecules.items() if not mols]
    if missing:
        raise RuntimeError(f"Queries without readable ligand SDFs: {missing}")

    rows = []
    target_read_failures = 0
    alignment_failures = 0
    score_failures = 0
    for target_id in target_ids:
        target_molecules, failures = read_ligands(
            args.target_systems_dir / target_id
        )
        target_read_failures += failures
        for query_id in query_ids:
            scores = []
            for query_mol in query_molecules[query_id]:
                for target_mol in target_molecules:
                    mobile = Chem.Mol(target_mol)
                    try:
                        sucos.align_molecules(query_mol, mobile)
                    except Exception:
                        alignment_failures += 1
                        continue
                    try:
                        scores.append(sucos.get_sucos_score(query_mol, mobile) * 100)
                    except Exception:
                        score_failures += 1
                        scores.append(0.0)
            rows.append(
                {
                    "query_system": query_id,
                    "target_system": target_id,
                    "sucos_shape": max(scores, default=np.nan),
                }
            )

    output = pd.DataFrame(rows)
    args.output_parquet.parent.mkdir(parents=True, exist_ok=True)
    output.to_parquet(args.output_parquet, index=False)
    print(f"Query systems: {len(query_ids):,}")
    print(f"Target systems: {len(target_ids):,}")
    print(f"Output rows: {len(output):,}")
    print(f"Calculated values: {output['sucos_shape'].notna().sum():,}")
    print(f"Query SDF read failures: {query_read_failures:,}")
    print(f"Target SDF read failures: {target_read_failures:,}")
    print(f"Alignment failures: {alignment_failures:,}")
    print(f"SuCOS failures recorded as zero: {score_failures:,}")
    print(f"Output: {args.output_parquet}")


if __name__ == "__main__":
    main()
