#!/bin/bash
#SBATCH --job-name=boltz2_reports
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/boltz2_report_%A_%a.out
#SBATCH --error=logs/boltz2_report_%A_%a.err

set -euo pipefail

PROJECT="${PLINDER_SIMILARITY_PROJECT:?Set PLINDER_SIMILARITY_PROJECT to the project directory}"
VERSION=boltz2_2023_06_01
ROOTFS="${PROJECT}/containers/plinder-rootfs"
APPTAINER="${APPTAINER_BIN:-apptainer}"
SCRIPT="${PROJECT}/scripts/similarity/consolidate_similarity_report.py"
MANIFEST="${PROJECT}/results/system_enumeration/calculable_queries/calculable_query_pdb_ids.txt"
PDB_ID="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST}" | tr -d '[:space:]' | tr '[:upper:]' '[:lower:]')"
[[ -n "${PDB_ID}" ]] || { echo "ERROR: empty PDB ID" >&2; exit 2; }

"${APPTAINER}" exec --cleanenv --env PATH=/opt/conda/bin:/usr/bin:/bin \
    --bind "${PROJECT}:/work" \
    --bind "${SCRIPT}:/report.py:ro" \
    "${ROOTFS}" /opt/conda/bin/python /report.py \
    --pdb-id "${PDB_ID}" \
    --plinder-scores "/work/results/${VERSION}/plinder/${PDB_ID}_combined_plinder_scores.parquet" \
    --morgan-scores "/work/results/${VERSION}/plinder/${PDB_ID}_morgan.parquet" \
    --sucos-scores "/work/results/${VERSION}/plinder/${PDB_ID}_sucos.parquet" \
    --reference-similarities /work/data/all_similarity_scores.parquet \
    --cutoff 2023-06-01 \
    --exclusive-cutoff \
    --output-dir "/work/results/${VERSION}/final_reports/${PDB_ID}"
