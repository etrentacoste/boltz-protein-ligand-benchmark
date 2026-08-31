#!/bin/bash
set -u

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 PDB_ID" >&2
    exit 2
fi

pdb_id="${1^^}"
root="$HOME/boltz_benchmark"
destination="$root/msa/${pdb_id}.csv"
generator="$root/scripts/generate_msa_only.sh"
max_attempts=5

if [[ -s "$destination" ]]; then
    echo "$pdb_id: validated destination MSA already exists"
    ls -lh "$destination"
    exit 0
fi

for attempt in $(seq 1 "$max_attempts"); do
    echo
    echo "========================================"
    echo "$pdb_id MSA attempt $attempt/$max_attempts"
    echo "========================================"

    temporary="$root/tmp/msa_generation_${pdb_id}"

    if [[ -d "$temporary" ]]; then
        failed="${temporary}_failed_$(date +%Y%m%d_%H%M%S)_pre${attempt}"
        mv "$temporary" "$failed"
        echo "Previous temporary directory saved: $failed"
    fi

    if [[ -e "$destination" ]]; then
        invalid="$root/msa/${pdb_id}_failed_$(date +%Y%m%d_%H%M%S)_attempt${attempt}.csv"
        mv "$destination" "$invalid"
        echo "Invalid destination moved to: $invalid"
    fi

    if "$generator" "$pdb_id"; then
        if [[ -s "$destination" ]]; then
            echo
            echo "$pdb_id MSA: SUCCESS on attempt $attempt"
            exit 0
        fi
    fi

    if [[ -d "$temporary" ]]; then
        failed="${temporary}_failed_$(date +%Y%m%d_%H%M%S)_attempt${attempt}"
        mv "$temporary" "$failed"
        echo "Failed attempt saved: $failed"
    fi

    if (( attempt < max_attempts )); then
        wait_seconds=$((attempt * 15))
        echo "Waiting ${wait_seconds}s before retry..."
        sleep "$wait_seconds"
    fi
done

echo "$pdb_id MSA: FAILED after $max_attempts attempts" >&2
exit 1
