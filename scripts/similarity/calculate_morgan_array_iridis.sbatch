#!/bin/bash
#SBATCH --job-name=plinder_morgan_array
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/morgan_array_%A_%a.out
#SBATCH --error=logs/morgan_array_%A_%a.err

set -euo pipefail

PROJECT="${PLINDER_SIMILARITY_PROJECT:?Set PLINDER_SIMILARITY_PROJECT to the project directory}"
ROOTFS="${PROJECT}/containers/plinder-rootfs"
APPTAINER="${APPTAINER_BIN:-apptainer}"
SCRIPT="${PROJECT}/scripts/similarity/calculate_morgan_array.py"
PDB_MANIFEST="${PROJECT}/results/system_enumeration/calculable_queries/calculable_query_pdb_ids.txt"
PDB_ID="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${PDB_MANIFEST}" | tr -d '[:space:]' | tr '[:upper:]' '[:lower:]')"
[[ -n "${PDB_ID}" ]] || { echo "ERROR: empty PDB ID" >&2; exit 2; }

"${APPTAINER}" exec --cleanenv \
    --env PATH=/opt/conda/bin:/usr/bin:/bin \
    --bind "${PROJECT}:/work" \
    --bind "${SCRIPT}:/calculate_morgan.py:ro" \
    "${ROOTFS}" /opt/conda/bin/python /calculate_morgan.py \
        --index-parquet /work/plinder/2024-06/v2/index/annotation_table.parquet \
        --query-systems-dir /work/queries/systems \
        --query-manifest "/work/results/query_curation/${PDB_ID}/selected_query_systems.txt" \
        --target-manifest "/work/results/system_enumeration/${PDB_ID}/calculable_target_systems.txt" \
        --output-parquet "/work/results/plinder/${PDB_ID}_morgan.parquet"
