#!/usr/bin/env python3

import argparse
import csv
import json
import shutil
import subprocess
import time
from pathlib import Path

import yaml


def validate_msa(path, expected_sequence):
    if not path.is_file() or path.stat().st_size < 20:
        return False, "missing_or_too_small"

    try:
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
    except Exception as error:
        return False, f"csv_error:{error}"

    if not rows:
        return False, "no_rows"

    if "sequence" not in rows[0]:
        return False, "missing_sequence_column"

    query = rows[0]["sequence"].replace("-", "").upper()

    if query != expected_sequence.upper():
        return False, (
            f"query_mismatch:{len(query)}_"
            f"vs_{len(expected_sequence)}"
        )

    return True, f"rows={len(rows)}"


def collect_sequences(batch_path, root):
    with batch_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    sequences = {}

    for row in rows:
        pdb_id = row["pdb_id"].upper()
        config_path = root / "config" / f"{pdb_id}.json"
        config = json.loads(config_path.read_text())

        for protein in config["proteins"]:
            sequence = protein["sequence"]
            digest = protein["sequence_hash"]

            if (
                digest in sequences
                and sequences[digest] != sequence
            ):
                raise RuntimeError(
                    f"Hash collision for {digest}"
                )

            sequences[digest] = sequence

    return sequences


def write_yaml(path, sequence):
    data = {
        "version": 1,
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": sequence,
                }
            }
        ],
    }

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        yaml.safe_dump(
            data,
            sort_keys=False,
        )
    )


def stop_process(process):
    if process.poll() is not None:
        return

    process.terminate()

    try:
        process.wait(timeout=15)
    except subprocess.TimeoutExpired:
        process.kill()
        process.wait(timeout=15)


def generate_one(
    digest,
    sequence,
    root,
    retries,
    timeout_minutes,
):
    destination = (
        root
        / "msa"
        / "by_sequence"
        / f"{digest}.csv"
    )
    destination.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    valid, detail = validate_msa(
        destination,
        sequence,
    )

    if valid:
        return "existing", detail

    yaml_path = (
        root
        / "inputs"
        / "msa_by_sequence"
        / f"{digest}.yaml"
    )
    write_yaml(yaml_path, sequence)

    log_dir = root / "logs" / "msa_by_sequence"
    log_dir.mkdir(parents=True, exist_ok=True)

    for attempt in range(1, retries + 1):
        attempt_dir = (
            root
            / "tmp"
            / "msa_by_sequence"
            / f"{digest}_attempt_{attempt}"
        )

        if attempt_dir.exists():
            shutil.rmtree(attempt_dir)

        attempt_dir.mkdir(parents=True)

        stdout_path = (
            log_dir
            / f"{digest}_attempt_{attempt}.out"
        )
        stderr_path = (
            log_dir
            / f"{digest}_attempt_{attempt}.err"
        )

        command = [
            "boltz",
            "predict",
            str(yaml_path),
            "--use_msa_server",
            "--accelerator",
            "cpu",
            "--devices",
            "1",
            "--recycling_steps",
            "1",
            "--sampling_steps",
            "1",
            "--diffusion_samples",
            "1",
            "--max_parallel_samples",
            "1",
            "--preprocessing-threads",
            "1",
            "--out_dir",
            str(attempt_dir),
        ]

        print(
            f"  attempt {attempt}/{retries}: "
            f"starting Boltz"
        )

        with stdout_path.open("w") as stdout_handle, \
             stderr_path.open("w") as stderr_handle:

            process = subprocess.Popen(
                command,
                stdout=stdout_handle,
                stderr=stderr_handle,
            )

            deadline = (
                time.time()
                + timeout_minutes * 60
            )
            last_candidate = None
            last_size = None
            stable_observations = 0
            copied = False
            copied_detail = ""

            while time.time() < deadline:
                candidates = sorted(
                    attempt_dir.rglob("*_0.csv")
                )

                for candidate in candidates:
                    valid, detail = validate_msa(
                        candidate,
                        sequence,
                    )

                    if not valid:
                        continue

                    size = candidate.stat().st_size

                    if (
                        candidate == last_candidate
                        and size == last_size
                    ):
                        stable_observations += 1
                    else:
                        last_candidate = candidate
                        last_size = size
                        stable_observations = 1

                    if stable_observations >= 3:
                        temporary = destination.with_suffix(
                            ".csv.part"
                        )
                        shutil.copy2(
                            candidate,
                            temporary,
                        )

                        valid_copy, copied_detail = (
                            validate_msa(
                                temporary,
                                sequence,
                            )
                        )

                        if not valid_copy:
                            temporary.unlink(
                                missing_ok=True
                            )
                            continue

                        temporary.replace(destination)
                        copied = True
                        break

                if copied:
                    stop_process(process)
                    break

                if process.poll() is not None:
                    break

                time.sleep(2)

            if not copied:
                stop_process(process)

        valid, final_detail = validate_msa(
            destination,
            sequence,
        )

        if valid:
            return "generated", final_detail

        print(
            f"  attempt {attempt} failed; "
            f"stdout={stdout_path} "
            f"stderr={stderr_path}"
        )

        time.sleep(3 * attempt)

    raise RuntimeError(
        f"Could not generate a valid MSA for {digest}"
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--batch",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--root",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--report",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=3,
    )
    parser.add_argument(
        "--timeout-minutes",
        type=int,
        default=20,
    )

    args = parser.parse_args()

    sequences = collect_sequences(
        args.batch,
        args.root,
    )

    print("Unique sequences:", len(sequences))

    results = []

    for index, (digest, sequence) in enumerate(
        sorted(sequences.items()),
        start=1,
    ):
        print()
        print(
            f"[{index}/{len(sequences)}] "
            f"{digest} length={len(sequence)}"
        )

        try:
            status, detail = generate_one(
                digest,
                sequence,
                args.root,
                args.retries,
                args.timeout_minutes,
            )

            result = {
                "sequence_hash": digest,
                "sequence_length": len(sequence),
                "status": status,
                "detail": detail,
                "msa_path": str(
                    args.root
                    / "msa"
                    / "by_sequence"
                    / f"{digest}.csv"
                ),
            }

            print(
                f"  SUCCESS: {status}, {detail}"
            )

        except Exception as error:
            result = {
                "sequence_hash": digest,
                "sequence_length": len(sequence),
                "status": "FAILED",
                "detail": str(error),
                "msa_path": "",
            }

            print(f"  FAILED: {error}")

        results.append(result)

    args.report.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with args.report.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sequence_hash",
                "sequence_length",
                "status",
                "detail",
                "msa_path",
            ],
        )
        writer.writeheader()
        writer.writerows(results)

    failures = [
        result
        for result in results
        if result["status"] == "FAILED"
    ]

    print()
    print("Sequences:", len(results))
    print("Successful:", len(results) - len(failures))
    print("Failed:", len(failures))
    print("Report:", args.report)

    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
