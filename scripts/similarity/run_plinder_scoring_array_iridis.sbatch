#!/bin/bash
#SBATCH --job-name=plinder_score_array
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/score_array_%A_%a.out
#SBATCH --error=logs/score_array_%A_%a.err

set -euo pipefail

PROJECT="${PLINDER_SIMILARITY_PROJECT:?Set PLINDER_SIMILARITY_PROJECT to the project directory}"
ROOTFS="${PROJECT}/containers/plinder-rootfs"
APPTAINER="${APPTAINER_BIN:-apptainer}"
SCRIPT="${PROJECT}/scripts/similarity/run_plinder_foldseek_scoring.py"
PDB_MANIFEST="${PROJECT}/results/system_enumeration/calculable_queries/calculable_query_pdb_ids.txt"

PDB_ID="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${PDB_MANIFEST}" | tr -d '[:space:]' | tr '[:upper:]' '[:lower:]')"
[[ -n "${PDB_ID}" ]] || { echo "ERROR: empty PDB ID" >&2; exit 2; }

QUERY_SYSTEMS="${PROJECT}/results/query_curation/${PDB_ID}/selected_query_systems.txt"
TARGET_SYSTEMS="${PROJECT}/results/system_enumeration/${PDB_ID}/calculable_target_systems.txt"
FOLDSEEK_TSV="${PROJECT}/results/foldseek/${PDB_ID}_foldseek_scoring.tsv"
MMSEQS_TSV="${PROJECT}/results/mmseqs/${PDB_ID}_mmseqs_scoring.tsv"
WORK_DIR="${PROJECT}/work/scoring_${PDB_ID}_${SLURM_JOB_ID}"
OUTPUT="${PROJECT}/results/plinder/${PDB_ID}_combined_plinder_scores.parquet"

for required in "${QUERY_SYSTEMS}" "${TARGET_SYSTEMS}" "${FOLDSEEK_TSV}" "${MMSEQS_TSV}"
do
    [[ -s "${required}" ]] || { echo "ERROR: required input missing: ${required}" >&2; exit 2; }
done

mkdir -p "${WORK_DIR}" "${PROJECT}/results/plinder"

"${APPTAINER}" exec \
    --cleanenv \
    --env "PLINDER_DIAGNOSE_MAPPING=${PLINDER_DIAGNOSE_MAPPING:-0}" \
    --env PATH=/opt/conda/bin:/usr/bin:/bin \
    --env PLINDER_RELEASE=2024-06 \
    --env PLINDER_ITERATION=v2 \
    --bind "${PROJECT}:/work" \
    --bind "${PROJECT}/plinder:/opt/plinder-data" \
    --bind "${SCRIPT}:/run_scoring.py:ro" \
    "${ROOTFS}" \
    /opt/conda/bin/python /run_scoring.py \
        --pdb-id "${PDB_ID}" \
        --query-system-manifest "/work/results/query_curation/${PDB_ID}/selected_query_systems.txt" \
        --foldseek-tsv "/work/results/foldseek/${PDB_ID}_foldseek_scoring.tsv" \
        --mmseqs-tsv "/work/results/mmseqs/${PDB_ID}_mmseqs_scoring.tsv" \
        --target-manifest "/work/results/system_enumeration/${PDB_ID}/calculable_target_systems.txt" \
        --plinder-dir /opt/plinder-data/2024-06/v2 \
        --entries-dir /work/plinder/2024-06/v2/reconstructed_entries \
        --work-dir "/work/work/scoring_${PDB_ID}_${SLURM_JOB_ID}" \
        --output-parquet "/work/results/plinder/${PDB_ID}_combined_plinder_scores.parquet"

if [[ -s "${OUTPUT}" ]]; then
    case "${WORK_DIR}" in
        "${PROJECT}/work/scoring_${PDB_ID}_"*) rm -rf -- "${WORK_DIR}" ;;
        *) echo "Refusing unsafe cleanup target: ${WORK_DIR}" >&2; exit 3 ;;
    esac
fi

echo "PDB ID: ${PDB_ID}"
echo "Output: ${OUTPUT}"
