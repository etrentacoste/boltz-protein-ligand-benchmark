#!/bin/bash
#SBATCH --job-name=boltz2_union
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/boltz2_union_%j.out
#SBATCH --error=logs/boltz2_union_%j.err

set -euo pipefail

PROJECT="${PLINDER_SIMILARITY_PROJECT:?Set PLINDER_SIMILARITY_PROJECT to the project directory}"
VERSION=boltz2_2023_06_01
ROOTFS="${PROJECT}/containers/plinder-rootfs"
APPTAINER="${APPTAINER_BIN:-apptainer}"
SCRIPT="${PROJECT}/scripts/similarity/summarize_multiquery_target_manifests.py"

"${APPTAINER}" exec --cleanenv \
    --bind "${PROJECT}:/work" \
    --bind "${SCRIPT}:/summarize.py:ro" \
    "${ROOTFS}" /opt/conda/bin/python /summarize.py \
    --query-manifest /work/results/system_enumeration/calculable_queries/calculable_query_pdb_ids.txt \
    --enumeration-dir "/work/results/${VERSION}/system_enumeration" \
    --systems-dir /work/plinder/2024-06/v2/systems \
    --output-dir "/work/results/${VERSION}/system_enumeration/all_queries"
