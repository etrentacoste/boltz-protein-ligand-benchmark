#!/usr/bin/env python3

import argparse
import csv
import gzip
import time
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def report_url(pdb_id):
    pdb_id = pdb_id.lower()
    middle = pdb_id[1:3]

    return (
        "https://files.rcsb.org/pub/pdb/"
        f"validation_reports/{middle}/{pdb_id}/"
        f"{pdb_id}_validation.xml.gz"
    )


def valid_report(path):
    if not path.is_file() or path.stat().st_size < 100:
        return False

    try:
        with gzip.open(path, "rb") as handle:
            ET.parse(handle)

        return True

    except Exception:
        return False


def download_one(pdb_id, output_dir, retries=3):
    pdb_id = pdb_id.upper()
    output = output_dir / f"{pdb_id}_validation.xml.gz"
    temporary = output_dir / f".{pdb_id}_validation.part"
    url = report_url(pdb_id)

    if valid_report(output):
        return pdb_id, "existing", output.stat().st_size

    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(
                url,
                headers={
                    "User-Agent": (
                        "boltz-benchmark-selection/1.0"
                    ),
                },
            )

            with urllib.request.urlopen(
                request,
                timeout=120,
            ) as response:
                content = response.read()

            temporary.write_bytes(content)

            if not valid_report(temporary):
                raise RuntimeError(
                    "Downloaded archive is not valid XML.gz"
                )

            temporary.replace(output)

            return (
                pdb_id,
                "downloaded",
                output.stat().st_size,
            )

        except urllib.error.HTTPError as error:
            if temporary.exists():
                temporary.unlink()

            if error.code == 404:
                return pdb_id, "not_available", 0

            if attempt == retries:
                return (
                    pdb_id,
                    f"failed_http_{error.code}",
                    0,
                )

            time.sleep(2 * attempt)

        except Exception as error:
            if temporary.exists():
                temporary.unlink()

            if attempt == retries:
                return pdb_id, f"failed:{error}", 0

            time.sleep(2 * attempt)

    return pdb_id, "failed_unknown", 0


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--report",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=6,
    )

    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    with args.input.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    pdb_ids = [
        row["pdb_id"].strip().upper()
        for row in rows
    ]

    results = {}

    with ThreadPoolExecutor(
        max_workers=args.workers
    ) as executor:
        futures = {
            executor.submit(
                download_one,
                pdb_id,
                args.output_dir,
            ): pdb_id
            for pdb_id in pdb_ids
        }

        completed = 0

        for future in as_completed(futures):
            pdb_id, status, size = future.result()
            results[pdb_id] = (status, size)
            completed += 1

            print(
                f"[{completed:3d}/{len(pdb_ids)}] "
                f"{pdb_id}: {status} "
                f"({size / 1024:.1f} KB)"
            )

    args.report.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with args.report.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "pdb_id",
            "status",
            "size_bytes",
        ])

        for pdb_id in pdb_ids:
            status, size = results[pdb_id]
            writer.writerow([
                pdb_id,
                status,
                size,
            ])

    available = sum(
        results[pdb_id][0] in {
            "downloaded",
            "existing",
        }
        for pdb_id in pdb_ids
    )

    not_available = sum(
        results[pdb_id][0] == "not_available"
        for pdb_id in pdb_ids
    )

    failed = (
        len(pdb_ids)
        - available
        - not_available
    )

    print()
    print("Requested:", len(pdb_ids))
    print("Available:", available)
    print("Not available:", not_available)
    print("Failed:", failed)
    print("Report:", args.report)

    if failed:
        print()
        print("Failed entries:")

        for pdb_id in pdb_ids:
            status, _ = results[pdb_id]

            if (
                status not in {
                    "downloaded",
                    "existing",
                    "not_available",
                }
            ):
                print(pdb_id, status)

        raise SystemExit(1)


if __name__ == "__main__":
    main()
