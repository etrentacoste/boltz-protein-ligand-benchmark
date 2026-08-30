#!/usr/bin/env python3

import argparse
import csv
import gzip
from pathlib import Path

import duckdb


REFERENCES = {
    "boltz1": {
        "cutoff": "2021-09-30",
        "expected_systems": 167_997,
        "expected_pdb_ids": 53_487,
    },
    "boltz2": {
        "cutoff": "2023-06-01",
        "expected_systems": 200_678,
        "expected_pdb_ids": 62_659,
    },
}


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Build the PLINDER-based Boltz-1 and Boltz-2 "
            "reference manifests used in the similarity analysis."
        )
    )

    parser.add_argument(
        "--similarity-parquet",
        required=True,
        type=Path,
        help=(
            "PLINDER-derived all_similarity_scores.parquet "
            "containing target_system and target_release_date."
        ),
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory in which manifests will be written.",
    )

    parser.add_argument(
        "--temp-dir",
        type=Path,
        default=None,
        help=(
            "Optional DuckDB temporary directory. Use a location "
            "with sufficient free space."
        ),
    )

    return parser.parse_args()


def systems_query(cutoff):
    return f"""
        SELECT DISTINCT
               target_system AS system_id
        FROM read_parquet(?)
        WHERE target_system IS NOT NULL
          AND trim(target_system) <> ''
          AND target_release_date < DATE '{cutoff}'
    """


def pdb_query(cutoff):
    reference_systems = systems_query(cutoff)

    return f"""
        SELECT DISTINCT
               upper(split_part(system_id, '__', 1)) AS pdb_id
        FROM ({reference_systems})
        WHERE split_part(system_id, '__', 1) <> ''
    """


def write_gzip_lines(
    connection,
    query,
    source_path,
    output_path,
):
    cursor = connection.execute(
        query,
        [str(source_path)],
    )

    number_written = 0

    with gzip.open(
        output_path,
        mode="wt",
        encoding="utf-8",
        newline="\n",
    ) as handle:
        while True:
            rows = cursor.fetchmany(10_000)

            if not rows:
                break

            for row in rows:
                handle.write(f"{row[0]}\n")

            number_written += len(rows)

    return number_written


def build_reference(
    connection,
    model_name,
    reference,
    source_path,
    output_dir,
):
    cutoff = reference["cutoff"]

    reference_systems_query = systems_query(cutoff)
    reference_pdb_query = pdb_query(cutoff)

    system_count = connection.execute(
        f"""
        SELECT count(*)
        FROM ({reference_systems_query})
        """,
        [str(source_path)],
    ).fetchone()[0]

    pdb_count = connection.execute(
        f"""
        SELECT count(*)
        FROM ({reference_pdb_query})
        """,
        [str(source_path)],
    ).fetchone()[0]

    if system_count != reference["expected_systems"]:
        raise RuntimeError(
            f"{model_name}: expected "
            f"{reference['expected_systems']} systems, "
            f"but found {system_count}."
        )

    if pdb_count != reference["expected_pdb_ids"]:
        raise RuntimeError(
            f"{model_name}: expected "
            f"{reference['expected_pdb_ids']} PDB IDs, "
            f"but found {pdb_count}."
        )

    minimum_date, maximum_date = connection.execute(
        f"""
        SELECT
            min(target_release_date),
            max(target_release_date)
        FROM read_parquet(?)
        WHERE target_system IS NOT NULL
          AND trim(target_system) <> ''
          AND target_release_date < DATE '{cutoff}'
        """,
        [str(source_path)],
    ).fetchone()

    systems_path = (
        output_dir
        / f"{model_name}_reference_systems.txt.gz"
    )

    pdb_ids_path = (
        output_dir
        / f"{model_name}_reference_pdb_ids.txt.gz"
    )

    written_systems = write_gzip_lines(
        connection,
        f"""
        SELECT system_id
        FROM ({reference_systems_query})
        ORDER BY system_id
        """,
        source_path,
        systems_path,
    )

    written_pdb_ids = write_gzip_lines(
        connection,
        f"""
        SELECT pdb_id
        FROM ({reference_pdb_query})
        ORDER BY pdb_id
        """,
        source_path,
        pdb_ids_path,
    )

    if written_systems != system_count:
        raise RuntimeError(
            f"{model_name}: wrote {written_systems} systems, "
            f"expected {system_count}."
        )

    if written_pdb_ids != pdb_count:
        raise RuntimeError(
            f"{model_name}: wrote {written_pdb_ids} PDB IDs, "
            f"expected {pdb_count}."
        )

    return {
        "model_name": model_name,
        "cutoff": cutoff,
        "system_count": system_count,
        "pdb_count": pdb_count,
        "minimum_date": minimum_date,
        "maximum_date": maximum_date,
        "systems_path": systems_path,
        "pdb_ids_path": pdb_ids_path,
    }


def write_summary(
    output_dir,
    source_path,
    result,
    extra_rows=None,
):
    model_name = result["model_name"]

    summary_path = (
        output_dir
        / f"{model_name}_reference_summary.csv"
    )

    rows = [
        (
            "reference_definition",
            f"PLINDER-based {model_name.capitalize()} "
            "reference universe",
        ),
        ("source_file", source_path.name),
        ("cutoff", result["cutoff"]),
        ("cutoff_operator", "<"),
        ("reference_systems", result["system_count"]),
        ("reference_pdb_ids", result["pdb_count"]),
    ]

    if extra_rows:
        rows.extend(extra_rows)

    rows.extend([
        (
            "minimum_release_date",
            result["minimum_date"],
        ),
        (
            "maximum_release_date",
            result["maximum_date"],
        ),
    ])

    with summary_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["item", "value"])
        writer.writerows(rows)

    return summary_path


def main():
    args = parse_arguments()

    source_path = args.similarity_parquet.resolve()
    output_dir = args.output_dir.resolve()

    if not source_path.is_file():
        raise FileNotFoundError(source_path)

    output_dir.mkdir(parents=True, exist_ok=True)

    if args.temp_dir is None:
        temporary_dir = output_dir / ".duckdb_tmp"
    else:
        temporary_dir = args.temp_dir.resolve()

    temporary_dir.mkdir(parents=True, exist_ok=True)

    connection = duckdb.connect()

    escaped_temp_directory = (
        str(temporary_dir).replace("'", "''")
    )

    connection.execute(
        f"SET temp_directory='{escaped_temp_directory}'"
    )

    results = {}

    for model_name, reference in REFERENCES.items():
        print()
        print("=" * 72)
        print(model_name)
        print("=" * 72)

        result = build_reference(
            connection=connection,
            model_name=model_name,
            reference=reference,
            source_path=source_path,
            output_dir=output_dir,
        )

        results[model_name] = result

        print("Cutoff:", f"< {result['cutoff']}")
        print("Systems:", result["system_count"])
        print("PDB IDs:", result["pdb_count"])

    boltz1_pdb_query = pdb_query(
        REFERENCES["boltz1"]["cutoff"]
    )

    boltz2_pdb_query = pdb_query(
        REFERENCES["boltz2"]["cutoff"]
    )

    overlap_count = connection.execute(
        f"""
        SELECT count(*)
        FROM ({boltz1_pdb_query}) AS boltz1
        INNER JOIN ({boltz2_pdb_query}) AS boltz2
        USING (pdb_id)
        """,
        [str(source_path), str(source_path)],
    ).fetchone()[0]

    boltz1_missing_from_boltz2 = (
        results["boltz1"]["pdb_count"]
        - overlap_count
    )

    additional_boltz2_pdb_ids = (
        results["boltz2"]["pdb_count"]
        - overlap_count
    )

    if boltz1_missing_from_boltz2 != 0:
        raise RuntimeError(
            f"{boltz1_missing_from_boltz2} Boltz-1 PDB IDs "
            "are absent from the Boltz-2 reference."
        )

    boltz1_summary = write_summary(
        output_dir=output_dir,
        source_path=source_path,
        result=results["boltz1"],
    )

    boltz2_summary = write_summary(
        output_dir=output_dir,
        source_path=source_path,
        result=results["boltz2"],
        extra_rows=[
            (
                "reused_boltz1_pdb_ids",
                overlap_count,
            ),
            (
                "additional_pdb_ids",
                additional_boltz2_pdb_ids,
            ),
            (
                "boltz1_ids_missing_from_boltz2",
                boltz1_missing_from_boltz2,
            ),
        ],
    )

    connection.close()

    print()
    print("=" * 72)
    print("REFERENCE MANIFEST GENERATION: SUCCESS")
    print("=" * 72)
    print("Boltz-1 summary:", boltz1_summary)
    print("Boltz-2 summary:", boltz2_summary)
    print("Output directory:", output_dir)


if __name__ == "__main__":
    main()
