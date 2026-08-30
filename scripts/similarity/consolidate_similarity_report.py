#!/usr/bin/env python3
"""Consolidate PLINDER, Morgan and SuCOS scores and create final reports."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


METRICS = (
    ("Protein sequence", "protein_fident_max", None),
    ("Protein structure", "protein_lddt_max", None),
    ("Pocket coverage", "pocket_qcov", 70.0),
    ("Protein-ligand interaction coverage", "pli_qcov", 50.0),
    ("Ligand pose (SuCOS shape)", "sucos_shape", 50.0),
    ("Ligand chemistry (Morgan Tanimoto)", "morgan_tanimoto", 50.0),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb-id", required=True)
    parser.add_argument("--plinder-scores", required=True, type=Path)
    parser.add_argument("--morgan-scores", required=True, type=Path)
    parser.add_argument("--sucos-scores", required=True, type=Path)
    parser.add_argument("--reference-similarities", required=True, type=Path)
    parser.add_argument("--cutoff", default="2021-09-30")
    parser.add_argument(
        "--exclusive-cutoff",
        action="store_true",
        help="Use release_date < cutoff instead of release_date <= cutoff.",
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def training_system_count(path: Path, cutoff: str, exclusive: bool) -> int:
    frame = pd.read_parquet(
        path, columns=["target_system", "target_release_date"]
    )
    dates = pd.to_datetime(frame["target_release_date"], errors="coerce")
    cutoff_date = pd.Timestamp(cutoff)
    selected = dates.lt(cutoff_date) if exclusive else dates.le(cutoff_date)
    return int(frame.loc[selected, "target_system"].nunique())


def load_plinder_wide(path: Path) -> pd.DataFrame:
    long = pd.read_parquet(
        path,
        columns=["query_system", "target_system", "metric", "similarity"],
    )
    wanted = {column for _, column, _ in METRICS[:4]}
    long = long[long["metric"].isin(wanted)]
    wide = (
        long.pivot_table(
            index=["query_system", "target_system"],
            columns="metric",
            values="similarity",
            aggfunc="max",
            observed=False,
        )
        .reset_index()
        .rename_axis(columns=None)
    )
    for column in wanted:
        if column not in wide:
            wide[column] = np.nan
    return wide


def format_value(value: float) -> str:
    return "not calculated" if pd.isna(value) else f"{value:.2f}%"


def main() -> None:
    args = parse_args()
    pdb_id = args.pdb_id.lower()
    total_training = training_system_count(
        args.reference_similarities, args.cutoff, args.exclusive_cutoff
    )

    detail = load_plinder_wide(args.plinder_scores)
    morgan = pd.read_parquet(args.morgan_scores)
    sucos = pd.read_parquet(args.sucos_scores)
    keys = ["query_system", "target_system"]
    detail = detail.merge(morgan, on=keys, how="outer")
    detail = detail.merge(sucos, on=keys, how="outer")
    detail.insert(0, "query_pdb_id", pdb_id)

    # A PDB-level report represents the complete user-supplied complex. When
    # curation retained multiple query systems, keep the best score per target.
    metric_columns = [column for _, column, _ in METRICS]
    per_target = (
        detail.groupby("target_system", as_index=False)[metric_columns]
        .max()
        .sort_values("target_system")
    )

    rows = []
    text_blocks = [f"PDB ID: {pdb_id.upper()}", f"Training systems: {total_training:,}"]
    for label, column, threshold in METRICS:
        values = pd.to_numeric(per_target[column], errors="coerce").dropna()
        calculated = int(len(values))
        calculated_pct = 100 * calculated / total_training if total_training else np.nan
        row = {
            "pdb_id": pdb_id,
            "metric": label,
            "metric_column": column,
            "threshold": threshold,
            "training_systems": total_training,
            "systems_calculated": calculated,
            "percent_training_calculated": calculated_pct,
            "minimum_calculated": values.min() if calculated else np.nan,
            "median_calculated": values.median() if calculated else np.nan,
            "maximum_calculated": values.max() if calculated else np.nan,
        }
        block = ["", label]
        if threshold is None:
            row.update(
                systems_above_threshold=np.nan,
                percent_calculated_above_threshold=np.nan,
                percent_training_above_threshold=np.nan,
                minimum_passing=np.nan,
                median_passing=np.nan,
                maximum_passing=np.nan,
            )
            block.extend(
                [
                    "Threshold: informative only",
                    f"Systems calculated: {calculated:,} ({calculated_pct:.3f}%)",
                    f"Minimum: {format_value(row['minimum_calculated'])}",
                    f"Median: {format_value(row['median_calculated'])}",
                    f"Maximum: {format_value(row['maximum_calculated'])}",
                ]
            )
        else:
            passing = values[values.ge(threshold)]
            passing_count = int(len(passing))
            pct_calculated = 100 * passing_count / calculated if calculated else np.nan
            pct_training = 100 * passing_count / total_training if total_training else np.nan
            row.update(
                systems_above_threshold=passing_count,
                percent_calculated_above_threshold=pct_calculated,
                percent_training_above_threshold=pct_training,
                minimum_passing=passing.min() if passing_count else np.nan,
                median_passing=passing.median() if passing_count else np.nan,
                maximum_passing=passing.max() if passing_count else np.nan,
            )
            block.extend(
                [
                    f"Threshold: {threshold:.0f}%",
                    f"Systems calculated: {calculated:,} ({calculated_pct:.3f}% of training systems)",
                    f"Minimum calculated: {format_value(row['minimum_calculated'])}",
                    f"Median calculated: {format_value(row['median_calculated'])}",
                    f"Maximum calculated: {format_value(row['maximum_calculated'])}",
                    f"Systems above threshold: {passing_count:,} ({pct_calculated:.3f}% of calculated systems; {pct_training:.3f}% of training systems)",
                    (
                        f"Minimum among passing systems: {format_value(row['minimum_passing'])}"
                        if passing_count
                        else "Passing-system statistics: no systems passed the threshold"
                    ),
                    *(
                        [
                            f"Median among passing systems: {format_value(row['median_passing'])}",
                            f"Maximum among passing systems: {format_value(row['maximum_passing'])}",
                        ]
                        if passing_count
                        else []
                    ),
                ]
            )
        rows.append(row)
        text_blocks.extend(block)

    summary = pd.DataFrame(rows)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    detail.to_parquet(args.output_dir / "pairwise_scores.parquet", index=False)
    per_target.to_parquet(args.output_dir / "pdb_target_scores.parquet", index=False)
    summary.to_csv(args.output_dir / "similarity_summary.csv", index=False)
    summary.to_parquet(args.output_dir / "similarity_summary.parquet", index=False)
    report = "\n".join(text_blocks) + "\n"
    (args.output_dir / "similarity_report.txt").write_text(report)
    print(report)
    print(f"Detailed pairs: {len(detail):,}")
    print(f"Unique targets: {len(per_target):,}")
    print(f"Output directory: {args.output_dir}")


if __name__ == "__main__":
    main()
