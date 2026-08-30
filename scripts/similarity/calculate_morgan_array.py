#!/usr/bin/env python3
"""Calculate the validated PLINDER/RnP Morgan similarity for one query PDB."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem


def read_manifest(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--index-parquet", required=True, type=Path)
    parser.add_argument("--query-systems-dir", required=True, type=Path)
    parser.add_argument("--query-manifest", required=True, type=Path)
    parser.add_argument("--target-manifest", required=True, type=Path)
    parser.add_argument("--output-parquet", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    query_ids = read_manifest(args.query_manifest)
    target_ids = read_manifest(args.target_manifest)
    generator = AllChem.GetMorganGenerator(
        radius=2, fpSize=2048, includeChirality=False
    )

    query_fps: dict[str, list] = defaultdict(list)
    query_failures = 0
    for system_id in query_ids:
        for path in sorted(
            (args.query_systems_dir / system_id / "ligand_files").glob("*.sdf")
        ):
            molecule = Chem.MolFromMolFile(str(path))
            if molecule is None:
                query_failures += 1
                continue
            query_fps[system_id].append(generator.GetFingerprint(molecule))

    missing_queries = [system_id for system_id in query_ids if not query_fps[system_id]]
    if missing_queries:
        raise RuntimeError(f"Queries without readable ligand SDFs: {missing_queries}")

    index = pd.read_parquet(
        args.index_parquet,
        columns=[
            "system_id",
            "ligand_instance_chain",
            "ligand_rdkit_canonical_smiles",
            "ligand_is_proper",
        ],
        filters=[
            ("system_id", "in", target_ids),
            ("ligand_is_proper", "==", True),
        ],
    )
    index = index.drop_duplicates(
        ["system_id", "ligand_instance_chain", "ligand_rdkit_canonical_smiles"]
    )

    target_fps: dict[str, list] = defaultdict(list)
    target_failures = 0
    for row in index.itertuples(index=False):
        smiles = row.ligand_rdkit_canonical_smiles
        if smiles is None or (isinstance(smiles, float) and np.isnan(smiles)):
            target_failures += 1
            continue
        molecule = Chem.MolFromSmiles(str(smiles))
        if molecule is None:
            target_failures += 1
            continue
        target_fps[str(row.system_id)].append(generator.GetFingerprint(molecule))

    rows = []
    for query_id in query_ids:
        for target_id in target_ids:
            scores = [
                DataStructs.TanimotoSimilarity(query_fp, target_fp) * 100
                for query_fp in query_fps[query_id]
                for target_fp in target_fps.get(target_id, [])
            ]
            rows.append(
                {
                    "query_system": query_id,
                    "target_system": target_id,
                    "morgan_tanimoto": max(scores, default=np.nan),
                }
            )

    output = pd.DataFrame(rows)
    args.output_parquet.parent.mkdir(parents=True, exist_ok=True)
    output.to_parquet(args.output_parquet, index=False)
    print(f"Query systems: {len(query_ids):,}")
    print(f"Target systems: {len(target_ids):,}")
    print(f"Output rows: {len(output):,}")
    print(f"Calculated values: {output['morgan_tanimoto'].notna().sum():,}")
    print(f"Query SDF failures: {query_failures:,}")
    print(f"Target SMILES failures: {target_failures:,}")
    print(f"Output: {args.output_parquet}")


if __name__ == "__main__":
    main()
