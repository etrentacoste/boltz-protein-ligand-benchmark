import argparse
import json
from pathlib import Path


HEADER = """#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={memory_gb}G
#SBATCH --time={time_limit}
#SBATCH --output={stdout}
#SBATCH --error={stderr}

set -euo pipefail

module purge
module load conda/python3
eval "$(conda shell.bash hook)"
conda activate {environment}

echo "Node: $(hostname)"
echo "Start: $(date)"
which python
which boltz
python --version

"""


def command(
    *,
    input_path,
    model,
    output_path,
    recycling_steps,
    sampling_steps,
    diffusion_samples,
    affinity_steps=None,
    affinity_samples=None,
):
    lines = [
        f"boltz predict {input_path} \\",
        f"    --model {model} \\",
        "    --accelerator cpu \\",
        "    --devices 1 \\",
        f"    --recycling_steps {recycling_steps} \\",
        f"    --sampling_steps {sampling_steps} \\",
        f"    --diffusion_samples {diffusion_samples} \\",
        "    --max_parallel_samples 1 \\",
    ]

    if affinity_steps is not None:
        lines.extend([
            f"    --sampling_steps_affinity {affinity_steps} \\",
            f"    --diffusion_samples_affinity {affinity_samples} \\",
        ])

    lines.extend([
        "    --seed 12345 \\",
        "    --override \\",
        f"    --out_dir {output_path}",
        "",
        'echo "Finish: $(date)"',
        "",
    ])

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=Path)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path.home() / "boltz_benchmark",
    )
    args = parser.parse_args()

    config = json.loads(args.config.read_text())
    pdb_id = config["pdb_id"].upper()
    root = args.root

    total_length = 0

    for protein in config["proteins"]:
        protein_ids = protein["id"]

        if isinstance(protein_ids, list):
            copy_count = len(protein_ids)
        else:
            copy_count = 1

        total_length += (
            len(protein["sequence"]) * copy_count
        )

    resources = config.get("resources", {})
    cpus = int(resources.get("cpus", 8))

    default_memory = 64 if total_length > 350 else 32
    memory_gb = int(
        resources.get("memory_gb", default_memory)
    )

    environment = (
        Path.home() / "envs/boltz2-conda"
    )

    jobs = [
        {
            "filename": f"pilot_{pdb_id}_boltz2.slurm",
            "job_name": f"pilot_b2_{pdb_id}",
            "model": "boltz2",
            "phase": "pilot",
            "time_limit": "02:00:00",
            "recycling_steps": 1,
            "sampling_steps": 50,
            "diffusion_samples": 1,
            "affinity_steps": 50,
            "affinity_samples": 1,
        },
        {
            "filename": f"production_{pdb_id}_boltz1.slurm",
            "job_name": f"prod_b1_{pdb_id}",
            "model": "boltz1",
            "phase": "production",
            "time_limit": "12:00:00",
            "recycling_steps": 3,
            "sampling_steps": 200,
            "diffusion_samples": 5,
            "affinity_steps": None,
            "affinity_samples": None,
        },
        {
            "filename": f"production_{pdb_id}_boltz2.slurm",
            "job_name": f"prod_b2_{pdb_id}",
            "model": "boltz2",
            "phase": "production",
            "time_limit": "12:00:00",
            "recycling_steps": 3,
            "sampling_steps": 200,
            "diffusion_samples": 5,
            "affinity_steps": 200,
            "affinity_samples": 5,
        },
    ]

    script_dir = root / "scripts/generated"
    script_dir.mkdir(parents=True, exist_ok=True)

    for job in jobs:
        phase = job["phase"]
        model = job["model"]

        log_dir = root / "logs" / phase
        result_dir = (
            root / "results" / phase / model / pdb_id
        )

        log_dir.mkdir(parents=True, exist_ok=True)
        result_dir.mkdir(parents=True, exist_ok=True)

        stdout = (
            log_dir
            / f"{model}_{pdb_id}_%j.out"
        )
        stderr = (
            log_dir
            / f"{model}_{pdb_id}_%j.err"
        )
        input_path = (
            root / f"inputs/{model}/{pdb_id}.yaml"
        )

        header = HEADER.format(
            job_name=job["job_name"],
            cpus=cpus,
            memory_gb=memory_gb,
            time_limit=job["time_limit"],
            stdout=stdout,
            stderr=stderr,
            environment=environment,
        )

        body = command(
            input_path=input_path,
            model=model,
            output_path=result_dir,
            recycling_steps=job["recycling_steps"],
            sampling_steps=job["sampling_steps"],
            diffusion_samples=job["diffusion_samples"],
            affinity_steps=job["affinity_steps"],
            affinity_samples=job["affinity_samples"],
        )

        output = script_dir / job["filename"]
        output.write_text(header + body)
        output.chmod(0o750)
        print("Written:", output)

    print("Protein residues:", total_length)
    print("CPUs:", cpus)
    print("Memory:", f"{memory_gb}G")


if __name__ == "__main__":
    main()
