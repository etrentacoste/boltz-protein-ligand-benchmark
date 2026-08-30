#!/bin/bash
#SBATCH --job-name=plinder_reports
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/report_array_%A_%a.out
#SBATCH --error=logs/report_array_%A_%a.err

set -euo pipefail

PROJECT="${PLINDER_SIMILARITY_PROJECT:?Set PLINDER_SIMILARITY_PROJECT to the project directory}"
ROOTFS="${PROJECT}/containers/plinder-rootfs"
APPTAINER="${APPTAINER_BIN:-apptainer}"
SCRIPT="${PROJECT}/scripts/similarity/consolidate_similarity_report.py"
PDB_MANIFEST="${PROJECT}/results/system_enumeration/calculable_queries/calculable_query_pdb_ids.txt"
PDB_ID="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${PDB_MANIFEST}" | tr -d '[:space:]' | tr '[:upper:]' '[:lower:]')"
[[ -n "${PDB_ID}" ]] || { echo "ERROR: empty PDB ID" >&2; exit 2; }

"${APPTAINER}" exec --cleanenv \
    --env PATH=/opt/conda/bin:/usr/bin:/bin \
    --bind "${PROJECT}:/work" \
    --bind "${SCRIPT}:/consolidate_report.py:ro" \
    "${ROOTFS}" /opt/conda/bin/python /consolidate_report.py \
        --pdb-id "${PDB_ID}" \
        --plinder-scores "/work/results/plinder/${PDB_ID}_combined_plinder_scores.parquet" \
        --morgan-scores "/work/results/plinder/${PDB_ID}_morgan.parquet" \
        --sucos-scores "/work/results/plinder/${PDB_ID}_sucos.parquet" \
        --reference-similarities /work/data/all_similarity_scores.parquet \
        --cutoff 2021-09-30 \
        --output-dir "/work/results/final_reports/${PDB_ID}"
