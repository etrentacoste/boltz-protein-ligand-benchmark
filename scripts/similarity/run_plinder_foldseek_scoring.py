#!/usr/bin/env python3
"""Run PLINDER pocket/PLI scoring from an existing Foldseek alignment TSV."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import pandas as pd
from plinder.data.utils.annotations.aggregate_annotations import Entry
from plinder.data.utils.annotations import get_similarity_scores as similarity_module
from plinder.data.utils.annotations.get_similarity_scores import Scorer


FOLDSEEK_COLUMNS = [
    "query",
    "target",
    "qlen",
    "fident",
    "alnlen",
    "qstart",
    "qend",
    "tstart",
    "tend",
    "evalue",
    "bits",
    "qcov",
    "tcov",
    "qaln",
    "taln",
    "lddt",
    "lddtfull",
]

MMSEQS_COLUMNS = [
    "query",
    "target",
    "qlen",
    "fident",
    "alnlen",
    "qstart",
    "qend",
    "tstart",
    "tend",
    "evalue",
    "bits",
    "qcov",
    "tcov",
    "qaln",
    "taln",
]


def normalize_identifier(identifier: str) -> str:
    """Normalize the PDB component while preserving the chain identifier."""
    text = str(identifier)
    pdb_id, separator, chain = text.partition("_")
    return f"{pdb_id.lower()}_{chain}" if separator else text.lower()


def normalize_query_identifier(
    identifier: str,
    pdb_id: str,
    default_chain: str | None,
) -> str:
    """Normalize query names even when they include an experimental_ prefix."""
    text = str(identifier).strip()
    if "_" in text:
        chain = text.rsplit("_", 1)[-1]
    elif default_chain is not None:
        chain = default_chain
    else:
        raise RuntimeError(
            f"Query identifier {text!r} has no chain suffix and the "
            f"selected PLINDER systems for {pdb_id} contain multiple chains"
        )
    return f"{pdb_id.lower()}_{chain}"


def pdb_from_identifier(identifier: str) -> str:
    return str(identifier).split("_", 1)[0].lower()


def requested_target_pdb_ids(target_manifest: Path) -> set[str]:
    return {
        line.strip().split("__", 1)[0].lower()
        for line in target_manifest.read_text().splitlines()
        if line.strip()
    }


def selected_query_author_chain(args: argparse.Namespace) -> str | None:
    """Return the sole selected protein author-chain ID, if unambiguous."""
    pdb_id = args.pdb_id.lower()
    entry = Entry.from_json(
        args.entries_dir / f"{pdb_id}.json",
        clear_non_pocket_residues=True,
        load_for_scoring=True,
    )
    selected = requested_query_systems(args)
    author_chains = {
        entry.chains[instance_chain.split(".", 1)[1]].auth_id
        for system_id in selected
        if system_id in entry.systems
        for instance_chain in entry.systems[system_id].protein_chains_asym_id
    }
    if len(author_chains) == 1:
        return next(iter(author_chains))
    return None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb-id", default="8flz")
    parser.add_argument("--query-system", action="append", dest="query_systems")
    parser.add_argument("--query-system-manifest", type=Path)
    parser.add_argument("--foldseek-tsv", required=True, type=Path)
    parser.add_argument("--mmseqs-tsv", type=Path)
    parser.add_argument("--target-manifest", required=True, type=Path)
    parser.add_argument("--plinder-dir", required=True, type=Path)
    parser.add_argument("--entries-dir", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--output-parquet", required=True, type=Path)
    return parser.parse_args()


def requested_query_systems(args: argparse.Namespace) -> set[str]:
    systems = set(args.query_systems or [])
    if args.query_system_manifest is not None:
        systems.update(
            line.strip()
            for line in args.query_system_manifest.read_text().splitlines()
            if line.strip()
        )
    if not systems:
        raise ValueError("At least one query system must be supplied")
    return systems


def install_reconstructed_entry_loader(entries_dir: Path) -> None:
    """Make PLINDER load locally reconstructed Entry JSON files."""
    def load_entries_from_json(
        data_dir: Path,
        pdb_ids: set[str],
        load_for_scoring: bool = True,
        **_: object,
    ) -> dict[str, Entry]:
        entries: dict[str, Entry] = {}
        missing = []
        for pdb_id in sorted({str(value).lower() for value in pdb_ids}):
            path = entries_dir / f"{pdb_id}.json"
            if not path.is_file():
                missing.append(pdb_id)
                continue
            entries[pdb_id] = Entry.from_json(
                path,
                clear_non_pocket_residues=True,
                load_for_scoring=load_for_scoring,
            )
        if missing:
            print(
                f"Entry cache misses: {len(missing):,} "
                f"(alignment rows for these PDB IDs will be ignored)"
            )
        print(f"Reconstructed entries loaded: {len(entries):,}")
        return entries

    similarity_module.load_entries_from_zips = load_entries_from_json


def prepare_alignment(args: argparse.Namespace) -> Path:
    alignment = pd.read_csv(
        args.foldseek_tsv,
        sep="\t",
        names=FOLDSEEK_COLUMNS,
        header=None,
    )
    if alignment.empty:
        raise RuntimeError(f"Foldseek table is empty: {args.foldseek_tsv}")

    pdb_id = args.pdb_id.lower()
    default_chain = selected_query_author_chain(args)
    alignment["query"] = alignment["query"].map(
        lambda value: normalize_query_identifier(value, pdb_id, default_chain)
    )
    alignment["target"] = alignment["target"].map(normalize_identifier)
    alignment["query_pdb_id"] = alignment["query"].map(pdb_from_identifier)
    alignment["target_pdb_id"] = alignment["target"].map(pdb_from_identifier)

    alignment = alignment[alignment["query_pdb_id"] == pdb_id].copy()
    target_pdb_ids = requested_target_pdb_ids(args.target_manifest)
    alignment = alignment[
        alignment["target_pdb_id"].isin(target_pdb_ids)
    ].copy()
    if alignment.empty:
        raise RuntimeError(
            f"No Foldseek rows remained for query PDB {pdb_id} "
            "after target filtering"
        )

    output = args.work_dir / "holo_foldseek" / "aln" / f"{pdb_id}.parquet"
    output.parent.mkdir(parents=True, exist_ok=True)
    alignment.to_parquet(output, index=False)

    print(f"Prepared Foldseek rows: {len(alignment):,}")
    print(f"Foldseek target PDB IDs: {alignment['target_pdb_id'].nunique():,}")
    print(f"Alignment Parquet: {output}")
    return output


def prepare_mmseqs_alignment(args: argparse.Namespace) -> Path | None:
    if args.mmseqs_tsv is None:
        return None

    alignment = pd.read_csv(args.mmseqs_tsv, sep="\t", header=0)
    missing_columns = sorted(set(MMSEQS_COLUMNS) - set(alignment.columns))
    if missing_columns:
        raise RuntimeError(f"MMseqs table is missing columns: {missing_columns}")
    alignment = alignment[MMSEQS_COLUMNS].copy()

    if alignment.empty:
        print(
            f"MMseqs produced no alignments for {args.pdb_id.lower()}; "
            "continuing with Foldseek-only mappings"
        )
        return None

    pdb_id = args.pdb_id.lower()
    default_chain = selected_query_author_chain(args)
    alignment["query"] = alignment["query"].map(
        lambda value: normalize_query_identifier(value, pdb_id, default_chain)
    )
    alignment["target"] = alignment["target"].map(normalize_identifier)
    alignment["query_pdb_id"] = alignment["query"].map(pdb_from_identifier)
    alignment["target_pdb_id"] = alignment["target"].map(pdb_from_identifier)

    alignment = alignment[alignment["query_pdb_id"] == pdb_id].copy()
    target_pdb_ids = requested_target_pdb_ids(args.target_manifest)
    alignment = alignment[
        alignment["target_pdb_id"].isin(target_pdb_ids)
    ].copy()
    if alignment.empty:
        print(
            f"No usable MMseqs rows remained for {pdb_id}; "
            "continuing with Foldseek-only mappings"
        )
        return None

    output = args.work_dir / "holo_mmseqs" / "aln" / f"{pdb_id}.parquet"
    output.parent.mkdir(parents=True, exist_ok=True)
    alignment.to_parquet(output, index=False)

    print(f"Prepared MMseqs rows: {len(alignment):,}")
    print(f"MMseqs target PDB IDs: {alignment['target_pdb_id'].nunique():,}")
    print(f"MMseqs alignment Parquet: {output}")
    return output


def main() -> None:
    args = parse_args()
    query_systems = requested_query_systems(args)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    prepare_alignment(args)
    prepare_mmseqs_alignment(args)
    install_reconstructed_entry_loader(args.entries_dir)

    scorer = Scorer(
        entries={},
        source_to_full_db_file={},
        db_dir=args.work_dir,
        scores_dir=args.work_dir / "scores",
        minimum_threshold=0,
        minimum_thresholds={
            "pli_qcov": 0,
            "pocket_qcov": 0,
            "pli_unique_qcov": 0,
        },
    )

    if os.environ.get("PLINDER_DIAGNOSE_MAPPING") == "1":
        alignment_path = (
            args.work_dir / "holo_foldseek" / "aln" / f"{args.pdb_id.lower()}.parquet"
        )
        alignment = pd.read_parquet(alignment_path)
        entry_ids = {args.pdb_id.lower()} | set(alignment["target_pdb_id"])
        scorer.entries = similarity_module.load_entries_from_zips(
            data_dir=args.plinder_dir,
            pdb_ids=entry_ids,
            load_for_scoring=True,
        )
        print("Running uncaught Foldseek mapping diagnostic", flush=True)
        scorer.map_alignment_df(
            alignment_path,
            aln_type="foldseek",
            search_db="holo",
        )
        raise RuntimeError("Diagnostic unexpectedly completed without an error")

    score_path = scorer.get_score_df(
        data_dir=args.plinder_dir,
        pdb_id=args.pdb_id.lower(),
        search_db="holo",
        overwrite=True,
    )
    if not score_path.is_file():
        raise RuntimeError(f"PLINDER did not create the score file: {score_path}")

    scores = pd.read_parquet(score_path)
    targets = {
        line.strip()
        for line in args.target_manifest.read_text().splitlines()
        if line.strip()
    }

    filtered = scores[
        (scores["query_system"].isin(query_systems))
        & (scores["target_system"].isin(targets))
    ].copy()

    args.output_parquet.parent.mkdir(parents=True, exist_ok=True)
    filtered.to_parquet(args.output_parquet, index=False)

    print(f"Raw PLINDER score rows: {len(scores):,}")
    print(f"Requested query systems: {len(query_systems):,}")
    print(f"Query systems represented: {filtered['query_system'].nunique():,}")
    print(f"Requested target systems: {len(targets):,}")
    print(f"Filtered query-target rows: {len(filtered):,}")
    print(f"Unique filtered targets: {filtered['target_system'].nunique():,}")

    report_columns = [
        column
        for column in (
            "protein_lddt_max",
            "protein_lddt_qcov_max",
            "pocket_qcov",
            "pli_qcov",
        )
        if column in filtered.columns
    ]
    print("Output metric columns:", report_columns)
    if report_columns:
        print(filtered[report_columns].describe().to_string())
    print(f"Output: {args.output_parquet}")


if __name__ == "__main__":
    main()
