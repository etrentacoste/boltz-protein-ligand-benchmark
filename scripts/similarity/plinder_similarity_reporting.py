"""Summarize query-to-training similarities for a PLINDER-like workflow.

This module does not calculate structural similarities. It consumes the Parquet
file produced by the scoring stage and creates an auditable report containing:

* the maximum score for each metric;
* the number of unique training systems at or above the metric threshold;
* the percentage among valid/evaluated training systems;
* the percentage relative to the complete training-set universe;
* the number of values that could not be calculated.

The same reporting code can be used first with Runs N' Poses scores for
validation and later with scores calculated for a new user-supplied PDB ID.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd


DEFAULT_CUTOFF = pd.Timestamp("2021-09-30")


@dataclass(frozen=True)
class MetricSpec:
    """Definition of one similarity metric and its optional threshold."""

    label: str
    column: str
    threshold: float | None


METRICS = (
    MetricSpec(
        label="Protein sequence",
        column="protein_fident_max",
        threshold=None,
    ),
    MetricSpec(
        label="Protein structure",
        column="protein_lddt_max",
        threshold=None,
    ),
    MetricSpec(
        label="Pocket",
        column="pocket_qcov",
        threshold=70.0,
    ),
    MetricSpec(
        label="Interaction",
        column="pli_qcov",
        threshold=50.0,
    ),
    MetricSpec(
        label="Ligand pose",
        column="sucos_shape",
        threshold=50.0,
    ),
    MetricSpec(
        label="Ligand chemistry",
        column="morgan_tanimoto",
        threshold=50.0,
    ),
)


def _normalise_system_ids(values: Iterable[object]) -> pd.Series:
    """Return clean, non-empty PLINDER system IDs."""

    series = pd.Series(values, dtype="string").dropna().str.strip()
    return series[series.ne("")]


def load_training_system_ids(
    similarity_file: Path,
    cutoff: pd.Timestamp = DEFAULT_CUTOFF,
) -> set[str]:
    """Extract the unique pre-cutoff target systems from Runs N' Poses.

    Notes
    -----
    ``all_similarity_scores.parquet`` contains targets released after the
    September 2021 training cutoff. Therefore, taking every unique
    ``target_system`` without applying ``target_release_date`` would include
    post-training structures.
    """

    similarity_file = Path(similarity_file)
    if not similarity_file.exists():
        raise FileNotFoundError(f"Similarity file not found: {similarity_file}")

    targets = pd.read_parquet(
        similarity_file,
        columns=["target_system", "target_release_date"],
    )
    targets["target_release_date"] = pd.to_datetime(
        targets["target_release_date"],
        errors="coerce",
    )
    targets = targets.loc[
        targets["target_release_date"].le(cutoff),
        "target_system",
    ]
    return set(_normalise_system_ids(targets))


def read_query_scores(
    score_file: Path,
    query_systems: Iterable[str] | None = None,
) -> pd.DataFrame:
    """Read only the columns required by the reporting stage.

    ``score_file`` can be the Runs N' Poses file for validation or a new
    Parquet produced by the custom scoring pipeline.
    """

    score_file = Path(score_file)
    if not score_file.exists():
        raise FileNotFoundError(f"Score file not found: {score_file}")

    available_columns = set(
        pd.read_parquet(score_file, columns=[]).columns
    )
    # Some parquet engines return no schema columns for columns=[].
    if not available_columns:
        import pyarrow.parquet as pq

        available_columns = set(pq.ParquetFile(score_file).schema_arrow.names)

    required = {
        "query_system",
        "target_system",
        *(metric.column for metric in METRICS),
    }
    missing = required - available_columns
    if missing:
        raise KeyError(
            "The score file is missing required columns: "
            f"{sorted(missing)}"
        )

    columns = sorted(required)
    for optional in ("group_key", "target_release_date"):
        if optional in available_columns:
            columns.append(optional)

    filters = None
    if query_systems is not None:
        clean_queries = sorted(set(_normalise_system_ids(query_systems)))
        if not clean_queries:
            raise ValueError("The query-system list is empty.")
        filters = [("query_system", "in", clean_queries)]

    scores = pd.read_parquet(
        score_file,
        columns=columns,
        filters=filters,
    )
    if scores.empty:
        raise ValueError("No score rows matched the requested query systems.")

    return scores


def filter_training_targets(
    scores: pd.DataFrame,
    training_system_ids: set[str],
) -> pd.DataFrame:
    """Retain only rows whose target belongs to the training-set universe."""

    filtered = scores.loc[
        scores["target_system"].isin(training_system_ids)
    ].copy()
    if filtered.empty:
        raise ValueError(
            "None of the scored targets belongs to the training-set universe."
        )
    return filtered


def collapse_to_unique_target_systems(
    scores: pd.DataFrame,
) -> pd.DataFrame:
    """Collapse duplicate ligand-pair rows to one row per training system.

    Runs N' Poses can contain more than one row for a query/target system pair.
    Counting rows would therefore overestimate the number of similar training
    systems. The maximum score per metric is retained for each unique target.
    """

    group_columns = ["query_system"]
    if "group_key" in scores.columns:
        group_columns.append("group_key")
    group_columns.append("target_system")

    metric_columns = [metric.column for metric in METRICS]
    return (
        scores.groupby(group_columns, dropna=False)[metric_columns]
        .max()
        .reset_index()
    )


def summarize_similar_training_systems(
    target_scores: pd.DataFrame,
    total_training_systems: int,
) -> pd.DataFrame:
    """Create one summary row per query unit and metric."""

    if total_training_systems <= 0:
        raise ValueError("total_training_systems must be positive.")

    unit_columns = ["query_system"]
    if "group_key" in target_scores.columns:
        unit_columns.append("group_key")

    rows: list[dict[str, object]] = []
    grouped = target_scores.groupby(unit_columns, dropna=False)

    for unit_key, unit_df in grouped:
        if not isinstance(unit_key, tuple):
            unit_key = (unit_key,)
        unit_values = dict(zip(unit_columns, unit_key))

        evaluated_targets = int(unit_df["target_system"].nunique())

        for metric in METRICS:
            values = pd.to_numeric(unit_df[metric.column], errors="coerce")
            valid = values.dropna()
            valid_count = int(valid.shape[0])
            missing_count = evaluated_targets - valid_count
            maximum = float(valid.max()) if valid_count else float("nan")

            if metric.threshold is None:
                similar_count: int | pd._libs.missing.NAType = pd.NA
                percent_valid = float("nan")
                percent_total = float("nan")
                threshold_passed: bool | pd._libs.missing.NAType = pd.NA
                status = "informative_only"
            else:
                similar_count = int(valid.ge(metric.threshold).sum())
                percent_valid = (
                    100.0 * similar_count / valid_count
                    if valid_count
                    else float("nan")
                )
                percent_total = (
                    100.0 * similar_count / total_training_systems
                )
                threshold_passed = similar_count > 0
                status = (
                    "similar_training_systems_found"
                    if threshold_passed
                    else "no_similar_training_system_found"
                )

            rows.append(
                {
                    **unit_values,
                    "metric": metric.label,
                    "metric_column": metric.column,
                    "threshold": metric.threshold,
                    "maximum_score": maximum,
                    "similar_training_system_count": similar_count,
                    "threshold_passed": threshold_passed,
                    "valid_scored_target_count": valid_count,
                    "missing_score_target_count": missing_count,
                    "evaluated_target_count": evaluated_targets,
                    "total_training_system_count": total_training_systems,
                    "percent_of_valid_targets_above_threshold": percent_valid,
                    "percent_of_all_training_targets_above_threshold": (
                        percent_total
                    ),
                    "status": status,
                }
            )

    summary = pd.DataFrame(rows)
    metric_order = {
        metric.label: position
        for position, metric in enumerate(METRICS)
    }
    summary["_metric_order"] = summary["metric"].map(metric_order)
    summary = (
        summary.sort_values(unit_columns + ["_metric_order"])
        .drop(columns="_metric_order")
        .reset_index(drop=True)
    )
    return summary


def build_similarity_report(
    reference_similarity_file: Path,
    query_score_file: Path,
    query_systems: Iterable[str] | None = None,
    cutoff: pd.Timestamp = DEFAULT_CUTOFF,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build the summary and unique-target detail tables."""

    training_ids = load_training_system_ids(
        reference_similarity_file,
        cutoff=cutoff,
    )
    raw_scores = read_query_scores(
        query_score_file,
        query_systems=query_systems,
    )
    training_scores = filter_training_targets(raw_scores, training_ids)
    target_scores = collapse_to_unique_target_systems(training_scores)
    summary = summarize_similar_training_systems(
        target_scores,
        total_training_systems=len(training_ids),
    )
    return summary, target_scores


def save_report(
    summary: pd.DataFrame,
    target_scores: pd.DataFrame,
    output_dir: Path,
    save_target_detail: bool = False,
) -> None:
    """Save compact report outputs."""

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    summary.to_csv(
        output_dir / "similarity_summary.csv",
        index=False,
    )
    summary.to_parquet(
        output_dir / "similarity_summary.parquet",
        index=False,
    )

    if save_target_detail:
        target_scores.to_parquet(
            output_dir / "unique_target_scores.parquet",
            index=False,
        )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Count similar training systems and report maximum similarity."
        )
    )
    parser.add_argument(
        "--reference-similarities",
        type=Path,
        required=True,
        help="Runs N' Poses all_similarity_scores.parquet.",
    )
    parser.add_argument(
        "--query-scores",
        type=Path,
        required=True,
        help="Parquet containing newly calculated query-to-target scores.",
    )
    parser.add_argument(
        "--query-system",
        action="append",
        dest="query_systems",
        help="Optional query system to retain; repeat for multiple systems.",
    )
    parser.add_argument(
        "--cutoff",
        default=str(DEFAULT_CUTOFF.date()),
        help="Training cutoff in YYYY-MM-DD format.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--save-target-detail",
        action="store_true",
        help="Also save one row per unique query/target system.",
    )
    return parser.parse_args()


def main() -> None:
    """Run the reporting pipeline from the command line."""

    args = parse_args()
    summary, target_scores = build_similarity_report(
        reference_similarity_file=args.reference_similarities,
        query_score_file=args.query_scores,
        query_systems=args.query_systems,
        cutoff=pd.Timestamp(args.cutoff),
    )
    save_report(
        summary,
        target_scores,
        output_dir=args.output_dir,
        save_target_detail=args.save_target_detail,
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
