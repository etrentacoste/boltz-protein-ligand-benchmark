#!/usr/bin/env python3
"""Calculate PoseBench/RnP SuCOS shape scores and validate against RnP."""

from __future__ import annotations

import argparse
import os
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, rdShapeAlign, rdShapeHelpers
from rdkit.Chem.FeatMaps import FeatMaps
from tqdm import tqdm


FDEF = AllChem.BuildFeatureFactory(
    os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
)
FEAT_MAP_PARAMS = {
    family: FeatMaps.FeatMapParams()
    for family in FDEF.GetFeatureFamilies()
}
PHARMACOPHORE_FEATURES = (
    "Donor",
    "Acceptor",
    "NegIonizable",
    "PosIonizable",
    "ZnBinder",
    "Aromatic",
    "Hydrophobe",
    "LumpedHydrophobe",
)

QUERY_MOL: Chem.Mol | None = None
SYSTEMS_DIR: Path | None = None


def read_molecule(path: Path) -> Chem.Mol | None:
    molecule = Chem.MolFromMolFile(str(path))
    if molecule is not None:
        return molecule
    try:
        return Chem.MolFromMolFile(
            str(path), sanitize=False, strictParsing=False
        )
    except Exception:
        return None


def align_molecules_crippen(
    reference: Chem.Mol, mobile: Chem.Mol, iterations: int = 100
) -> None:
    alignment = Chem.rdMolAlign.GetCrippenO3A(
        mobile, reference, maxIters=iterations
    )
    alignment.Align()


def align_molecules(
    reference: Chem.Mol,
    mobile: Chem.Mol,
    max_preiters: int = 100,
    max_postiters: int = 100,
) -> tuple[float, float]:
    align_molecules_crippen(reference, mobile, iterations=max_preiters)
    return rdShapeAlign.AlignMol(
        reference,
        mobile,
        max_preiters=max_preiters,
        max_postiters=max_postiters,
    )


def get_feature_map_score(
    molecule_1: Chem.Mol,
    molecule_2: Chem.Mol,
    score_mode: FeatMaps.FeatMapScoreMode = FeatMaps.FeatMapScoreMode.All,
) -> float:
    feature_lists = []
    for molecule in (molecule_1, molecule_2):
        raw_features = FDEF.GetFeaturesForMol(molecule)
        feature_lists.append(
            [
                feature
                for feature in raw_features
                if feature.GetFamily() in PHARMACOPHORE_FEATURES
            ]
        )

    feature_maps = [
        FeatMaps.FeatMap(
            feats=features,
            weights=[1] * len(features),
            params=FEAT_MAP_PARAMS,
        )
        for features in feature_lists
    ]
    feature_maps[0].scoreMode = score_mode
    denominator = min(feature_maps[0].GetNumFeatures(), len(feature_lists[1]))
    if denominator == 0:
        raise ValueError("No comparable pharmacophore features")
    return feature_maps[0].ScoreFeats(feature_lists[1]) / denominator


def get_sucos_score(molecule_1: Chem.Mol, molecule_2: Chem.Mol) -> float:
    feature_score = np.clip(
        get_feature_map_score(molecule_1, molecule_2), 0, 1
    )
    protrude_distance = np.clip(
        rdShapeHelpers.ShapeProtrudeDist(
            molecule_1, molecule_2, allowReordering=False
        ),
        0,
        1,
    )
    return float(0.5 * feature_score + 0.5 * (1 - protrude_distance))


def initialize_worker(query_sdf: str, systems_dir: str) -> None:
    global QUERY_MOL, SYSTEMS_DIR
    QUERY_MOL = read_molecule(Path(query_sdf))
    SYSTEMS_DIR = Path(systems_dir)
    if QUERY_MOL is None:
        raise RuntimeError(f"Could not read query ligand: {query_sdf}")


def score_system(system_id: str) -> dict[str, str | int | float]:
    if QUERY_MOL is None or SYSTEMS_DIR is None:
        raise RuntimeError("Worker was not initialized")

    ligand_paths = sorted(
        (SYSTEMS_DIR / system_id / "ligand_files").glob("*.sdf")
    )
    scores = []
    read_failures = 0
    alignment_failures = 0
    sucos_failures = 0

    for ligand_path in ligand_paths:
        mobile = read_molecule(ligand_path)
        if mobile is None:
            read_failures += 1
            continue
        try:
            align_molecules(QUERY_MOL, mobile)
        except Exception:
            alignment_failures += 1
            continue
        try:
            scores.append(get_sucos_score(QUERY_MOL, mobile))
        except Exception:
            sucos_failures += 1
            scores.append(0.0)

    return {
        "target_system": system_id,
        "ligand_files": len(ligand_paths),
        "read_failures": read_failures,
        "alignment_failures": alignment_failures,
        "sucos_failures": sucos_failures,
        "sucos_shape_calculated": max(scores, default=np.nan) * 100,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--systems-dir", required=True, type=Path)
    parser.add_argument("--query-system", required=True)
    parser.add_argument("--query-ligand", required=True)
    parser.add_argument("--target-manifest", required=True, type=Path)
    parser.add_argument("--rnp-parquet", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--output-parquet", required=True, type=Path)
    parser.add_argument("--comparison-csv", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    query_sdf = (
        args.systems_dir
        / args.query_system
        / "ligand_files"
        / args.query_ligand
    )
    target_ids = [
        line.strip()
        for line in args.target_manifest.read_text().splitlines()
        if line.strip()
    ]

    with Pool(
        processes=args.workers,
        initializer=initialize_worker,
        initargs=(str(query_sdf), str(args.systems_dir)),
    ) as pool:
        rows = list(
            tqdm(
                pool.imap_unordered(score_system, target_ids, chunksize=10),
                total=len(target_ids),
                desc="SuCOS shape",
            )
        )

    calculated = pd.DataFrame(rows).sort_values("target_system")
    args.output_parquet.parent.mkdir(parents=True, exist_ok=True)
    calculated.to_parquet(args.output_parquet, index=False)

    rnp = pd.read_parquet(
        args.rnp_parquet,
        columns=["query_system", "target_system", "sucos_shape"],
    )
    rnp = rnp[rnp["query_system"] == args.query_system].drop_duplicates(
        "target_system"
    )
    comparison = calculated.merge(rnp, on="target_system", how="left")
    comparison.to_csv(args.comparison_csv, index=False)

    valid = comparison.dropna(
        subset=["sucos_shape_calculated", "sucos_shape"]
    )
    difference = (
        valid["sucos_shape_calculated"] - valid["sucos_shape"]
    ).abs()
    exact = np.isclose(
        valid["sucos_shape_calculated"], valid["sucos_shape"], atol=1e-8
    )
    close_01 = np.isclose(
        valid["sucos_shape_calculated"], valid["sucos_shape"], atol=0.1
    )

    print(f"Requested systems: {len(target_ids):,}")
    print(
        "Calculated systems:",
        f"{comparison['sucos_shape_calculated'].notna().sum():,}",
    )
    print(f"Comparable systems: {len(valid):,}")
    print(f"Read failures: {calculated['read_failures'].sum():,}")
    print(f"Alignment failures: {calculated['alignment_failures'].sum():,}")
    print(f"SuCOS failures recorded as zero: {calculated['sucos_failures'].sum():,}")
    if len(valid):
        print(f"Exact matches: {exact.sum():,} ({exact.mean() * 100:.2f}%)")
        print(
            f"Matches within 0.1: {close_01.sum():,} "
            f"({close_01.mean() * 100:.2f}%)"
        )
        print(
            "Pearson correlation:",
            f"{valid['sucos_shape_calculated'].corr(valid['sucos_shape']):.6f}",
        )
        print(f"Mean absolute difference: {difference.mean():.6f}")
        print(f"Median absolute difference: {difference.median():.6f}")
        print(f"Maximum absolute difference: {difference.max():.6f}")
        print(
            "Calculated >= 50:",
            f"{(valid['sucos_shape_calculated'] >= 50).sum():,}",
        )
        print(f"RnP >= 50: {(valid['sucos_shape'] >= 50).sum():,}")
    print(f"Output: {args.output_parquet}")
    print(f"Comparison: {args.comparison_csv}")


if __name__ == "__main__":
    main()

