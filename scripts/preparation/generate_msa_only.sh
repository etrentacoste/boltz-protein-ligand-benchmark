#!/bin/bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 PDB_ID" >&2
    exit 2
fi

pdb_id="${1^^}"
root="$HOME/boltz_benchmark"
input="$root/inputs/msa_generation/${pdb_id}.yaml"
temporary="$root/tmp/msa_generation_${pdb_id}"
destination="$root/msa/${pdb_id}.csv"
log="$root/logs/msa_${pdb_id}.log"

if [[ ! -s "$input" ]]; then
    echo "ERROR: missing input $input" >&2
    exit 1
fi

mkdir -p "$temporary" "$root/msa" "$root/logs"

if [[ -s "$destination" ]]; then
    echo "MSA already exists: $destination"
    exit 0
fi

echo "Starting MSA generation for $pdb_id"
echo "Log: $log"

# Start Boltz in a separate process group. The process will be stopped
# as soon as the final MSA CSV has been written.
setsid boltz predict "$input" \
    --use_msa_server \
    --accelerator cpu \
    --devices 1 \
    --recycling_steps 1 \
    --sampling_steps 50 \
    --diffusion_samples 1 \
    --max_parallel_samples 1 \
    --preprocessing-threads 1 \
    --override \
    --out_dir "$temporary" \
    >"$log" 2>&1 &

boltz_pid=$!
echo "Boltz PID: $boltz_pid"

msa_file=""
previous_size=0
stable_checks=0

for attempt in $(seq 1 900); do
    msa_file="$(
        find "$temporary" \
            -type f \
            -path "*/msa/*_0.csv" \
            -size +100c \
            -print -quit 2>/dev/null || true
    )"

    if [[ -n "$msa_file" ]]; then
        current_size=$(stat -c '%s' "$msa_file")

        if [[ "$current_size" -eq "$previous_size" ]]; then
            stable_checks=$((stable_checks + 1))
        else
            stable_checks=0
        fi

        previous_size="$current_size"

        if [[ "$stable_checks" -ge 2 ]]; then
            break
        fi
    fi

    if ! kill -0 "$boltz_pid" 2>/dev/null; then
        echo "ERROR: Boltz stopped before a stable MSA was found." >&2
        tail -n 50 "$log" >&2
        exit 1
    fi

    sleep 2
done

if [[ -z "$msa_file" || "$stable_checks" -lt 2 ]]; then
    echo "ERROR: MSA generation timed out." >&2
    kill -TERM -- "-$boltz_pid" 2>/dev/null || true
    wait "$boltz_pid" 2>/dev/null || true
    exit 1
fi

cp "$msa_file" "$destination"

echo "Final MSA detected: $msa_file"
echo "Copied to: $destination"

# Stop the entire Boltz process group before structural inference consumes
# login-node resources.
kill -TERM -- "-$boltz_pid" 2>/dev/null || true
wait "$boltz_pid" 2>/dev/null || true

python - "$pdb_id" "$destination" "$root/config/${pdb_id}.json" <<'PY_VALIDATE'
import csv
import json
import sys
from pathlib import Path

pdb_id = sys.argv[1]
msa_path = Path(sys.argv[2])
config_path = Path(sys.argv[3])

config = json.loads(config_path.read_text())
expected = config["proteins"][0]["sequence"]

with msa_path.open() as handle:
    rows = list(csv.DictReader(handle))

if not rows:
    raise SystemExit("ERROR: MSA contains no sequences")

query = rows[0]["sequence"].replace("-", "")

print("Columns:", rows[0].keys())
print("MSA rows:", len(rows))
print("Query length:", len(query))
print("Expected length:", len(expected))
print("Query matches input:", query == expected)

if query != expected:
    raise SystemExit(
        f"ERROR: MSA query does not match {pdb_id}"
    )

print("MSA validation: SUCCESS")
PY_VALIDATE

echo "MSA-only workflow completed for $pdb_id"
