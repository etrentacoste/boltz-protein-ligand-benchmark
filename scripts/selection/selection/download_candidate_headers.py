#!/usr/bin/env python3

import argparse
import csv
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def download_one(pdb_id, output_dir, retries=3):
    pdb_id = pdb_id.upper()
    output = output_dir / f"{pdb_id}.cif"
    url = f"https://files.rcsb.org/header/{pdb_id}.cif"

    if output.is_file() and output.stat().st_size > 100:
        text = output.read_text(errors="replace")

        if text.lstrip().startswith("data_"):
            return pdb_id, "existing", output.stat().st_size

    temporary = output.with_suffix(".cif.part")

    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(
                url,
                headers={
                    "User-Agent": "boltz-benchmark-selection/1.0",
                },
            )

            with urllib.request.urlopen(
                request,
                timeout=60,
            ) as response:
                content = response.read()

            text = content.decode("utf-8", errors="replace")

            if len(content) < 100:
                raise RuntimeError(
                    f"Downloaded file is too small: {len(content)} bytes"
                )

            if not text.lstrip().startswith("data_"):
                raise RuntimeError(
                    "Downloaded content is not an mmCIF file"
                )

            temporary.write_bytes(content)
            temporary.replace(output)

            return pdb_id, "downloaded", len(content)

        except Exception as error:
            if temporary.exists():
                temporary.unlink()

            if attempt == retries:
                return pdb_id, f"failed: {error}", 0

            time.sleep(2 * attempt)

    return pdb_id, "failed: unknown error", 0


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
        "--workers",
        type=int,
        default=8,
    )

    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    with args.input.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    pdb_ids = [
        row["pdb_id"].strip().upper()
        for row in rows
    ]

    results = []

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
            results.append((pdb_id, status, size))
            completed += 1

            print(
                f"[{completed:3d}/{len(pdb_ids)}] "
                f"{pdb_id}: {status} ({size} bytes)"
            )

    result_by_id = {
        pdb_id: (status, size)
        for pdb_id, status, size in results
    }

    report = args.output_dir.parent / "header_download_report.csv"

    with report.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["pdb_id", "status", "size_bytes"])

        for pdb_id in pdb_ids:
            status, size = result_by_id[pdb_id]
            writer.writerow([pdb_id, status, size])

    failures = [
        (pdb_id, status)
        for pdb_id, status, _ in results
        if status.startswith("failed")
    ]

    print()
    print("Requested:", len(pdb_ids))
    print("Successful:", len(pdb_ids) - len(failures))
    print("Failed:", len(failures))
    print("Report:", report)

    if failures:
        print()
        print("Failed identifiers:")

        for pdb_id, status in sorted(failures):
            print(" ", pdb_id, status)

        raise SystemExit(1)


if __name__ == "__main__":
    main()
