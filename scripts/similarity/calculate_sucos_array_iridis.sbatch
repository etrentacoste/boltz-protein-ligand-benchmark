#!/bin/bash
#SBATCH --job-name=plinder_sucos_array
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/sucos_array_%A_%a.out
#SBATCH --error=logs/sucos_array_%A_%a.err

set -euo pipefail

PROJECT="${PLINDER_SIMILARITY_PROJECT:?Set PLINDER_SIMILARITY_PROJECT to the project directory}"
PYTHON="${PROJECT}/envs/sucos-rdkit2024/bin/python"
SCRIPT_DIR="${PROJECT}/scripts/similarity"
PDB_MANIFEST="${PROJECT}/results/system_enumeration/calculable_queries/calculable_query_pdb_ids.txt"
PDB_ID="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${PDB_MANIFEST}" | tr -d '[:space:]' | tr '[:upper:]' '[:lower:]')"
[[ -n "${PDB_ID}" ]] || { echo "ERROR: empty PDB ID" >&2; exit 2; }

export PYTHONPATH="${SCRIPT_DIR}:${PYTHONPATH:-}"
"${PYTHON}" "${SCRIPT_DIR}/calculate_sucos_array.py" \
    --query-systems-dir "${PROJECT}/queries/systems" \
    --target-systems-dir "${PROJECT}/plinder/2024-06/v2/systems" \
    --query-manifest "${PROJECT}/results/query_curation/${PDB_ID}/selected_query_systems.txt" \
    --target-manifest "${PROJECT}/results/system_enumeration/${PDB_ID}/calculable_target_systems.txt" \
    --output-parquet "${PROJECT}/results/plinder/${PDB_ID}_sucos.parquet"
