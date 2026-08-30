#!/usr/bin/env python3

import argparse
import csv
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def is_valid_cif(path):
    if not path.is_file() or path.stat().st_size < 1000:
        return False

    text = path.read_text(
        errors="replace",
    )

    return (
        text.lstrip().startswith("data_")
        and "_atom_site.Cartn_x" in text
        and "_atom_site.Cartn_y" in text
        and "_atom_site.Cartn_z" in text
    )


def download_one(pdb_id, output_dir, retries=3):
    pdb_id = pdb_id.upper()
    output = output_dir / f"{pdb_id}.cif"
    temporary = output.with_suffix(".cif.part")
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"

    if is_valid_cif(output):
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

            if not is_valid_cif(temporary):
                raise RuntimeError(
                    "Downloaded file is not a complete mmCIF"
                )

            temporary.replace(output)

            return (
                pdb_id,
                "downloaded",
                output.stat().st_size,
            )

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
        default=6,
    )
    parser.add_argument(
        "--report",
        type=Path,
        required=True,
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
            writer.writerow([pdb_id, status, size])

    failures = [
        pdb_id
        for pdb_id in pdb_ids
        if results[pdb_id][0].startswith("failed")
    ]

    print()
    print("Requested:", len(pdb_ids))
    print("Successful:", len(pdb_ids) - len(failures))
    print("Failed:", len(failures))
    print("Report:", args.report)

    if failures:
        print("Failed PDB IDs:", " ".join(failures))
        raise SystemExit(1)


if __name__ == "__main__":
    main()
