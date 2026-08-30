#!/bin/bash
#SBATCH --job-name=prod_b2_31WA
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --chdir=/home/et3n25/boltz_benchmark
#SBATCH --output=/home/et3n25/boltz_benchmark/logs/production/b2_31WA_%j.out
#SBATCH --error=/home/et3n25/boltz_benchmark/logs/production/b2_31WA_%j.err

set -euo pipefail

module purge
module load conda/python3
eval "$(conda shell.bash hook)"
conda activate /home/et3n25/envs/boltz2-conda

echo "Node: $(hostname)"
echo "Start: $(date)"
echo "Job ID: $SLURM_JOB_ID"
which python
which boltz
python --version

boltz predict /home/et3n25/boltz_benchmark/inputs/boltz2/31WA.yaml \
    --model boltz2 \
    --accelerator cpu \
    --devices 1 \
    --recycling_steps 3 \
    --sampling_steps 200 \
    --diffusion_samples 5 \
    --max_parallel_samples 1 \
    --sampling_steps_affinity 200 \
    --diffusion_samples_affinity 5 \
    --seed 12345 \
    --override \
    --out_dir /home/et3n25/boltz_benchmark/results/production/boltz2/31WA

echo "Finish: $(date)"
