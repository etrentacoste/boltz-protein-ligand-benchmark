#!/usr/bin/env python3

import argparse
import csv
import json
import urllib.request
from pathlib import Path


SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"


def terminal(attribute, operator, value):
    return {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": attribute,
            "operator": operator,
            "negation": False,
            "value": value,
        },
    }


def build_query(date_from, date_to):
    nodes = [
        terminal(
            "rcsb_accession_info.initial_release_date",
            "range",
            {
                "from": date_from,
                "to": date_to,
                "include_lower": True,
                "include_upper": True,
            },
        ),
        terminal(
            "exptl.method",
            "exact_match",
            "X-RAY DIFFRACTION",
        ),
        terminal(
            "rcsb_entry_info.resolution_combined",
            "less_or_equal",
            2.5,
        ),
        terminal(
            "rcsb_entry_info.polymer_entity_count_protein",
            "greater_or_equal",
            1,
        ),
        terminal(
            "rcsb_entry_info.nonpolymer_entity_count",
            "greater_or_equal",
            1,
        ),
        terminal(
            "rcsb_entity_source_organism.ncbi_scientific_name",
            "exact_match",
            "Homo sapiens",
        ),
        terminal(
            "entity_poly.rcsb_entity_polymer_type",
            "exact_match",
            "Protein",
        ),
    ]

    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": nodes,
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": 10000,
            },
            "results_content_type": ["experimental"],
            "sort": [
                {
                    "sort_by": (
                        "rcsb_accession_info.initial_release_date"
                    ),
                    "direction": "desc",
                }
            ],
        },
    }


def submit(query):
    payload = json.dumps(query).encode("utf-8")

    request = urllib.request.Request(
        SEARCH_URL,
        data=payload,
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": "boltz-benchmark-selection/1.0",
        },
        method="POST",
    )

    with urllib.request.urlopen(request, timeout=120) as response:
        return json.load(response)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--date-from",
        default="2025-01-10",
    )
    parser.add_argument(
        "--date-to",
        required=True,
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=(
            Path.home()
            / "boltz_benchmark"
            / "selection"
        ),
    )

    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    query = build_query(args.date_from, args.date_to)

    query_path = args.output_dir / "stage0_query.json"
    query_path.write_text(
        json.dumps(query, indent=2) + "\n"
    )

    result = submit(query)

    identifiers = [
        item["identifier"].upper()
        for item in result.get("result_set", [])
    ]

    output_path = (
        args.output_dir
        / "stage0_candidate_ids.csv"
    )

    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["pdb_id"])

        for pdb_id in identifiers:
            writer.writerow([pdb_id])

    print("Date range:", args.date_from, "to", args.date_to)
    print("Reported total:", result.get("total_count"))
    print("Identifiers written:", len(identifiers))
    print("Query:", query_path)
    print("Candidates:", output_path)

    if len(identifiers) != len(set(identifiers)):
        raise RuntimeError("Duplicate PDB identifiers found")

    if len(identifiers) == 10000:
        print(
            "WARNING: the result reached the pagination limit."
        )


if __name__ == "__main__":
    main()
