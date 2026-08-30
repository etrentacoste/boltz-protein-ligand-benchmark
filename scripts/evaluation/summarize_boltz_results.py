import argparse
import csv
import json
import statistics
from pathlib import Path


def load_json(path):
    with path.open() as handle:
        return json.load(handle)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb-id", required=True)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path.home() / "boltz_visualization",
    )
    args = parser.parse_args()

    pdb_id = args.pdb_id.upper()
    complex_dir = args.root / pdb_id
    rows = []
    summaries = []

    for model in ["boltz1", "boltz2"]:
        model_dir = complex_dir / "production" / model
        metrics_dir = complex_dir / "metrics" / model
        model_rows = []

        for pose in range(5):
            metric_path = metrics_dir / f"model_{pose}.json"
            confidence_path = (
                model_dir
                / f"confidence_{pdb_id}_model_{pose}.json"
            )

            metrics = load_json(metric_path)

            rmsd_data = metrics["rmsd"]["assigned_scores"][0]
            lddt_data = metrics["lddt_pli"]["assigned_scores"][0]

            rmsd = float(rmsd_data["score"])
            lddt = float(lddt_data["score"])
            coverage = float(lddt_data.get("coverage", 1.0))

            confidence = (
                load_json(confidence_path)
                if confidence_path.is_file()
                else {}
            )

            rmsd_success = rmsd < 2.0
            lddt_success = lddt > 0.8
            combined_success = rmsd_success and lddt_success

            row = {
                "pdb_id": pdb_id,
                "model": model,
                "pose": pose,
                "is_top1": pose == 0,
                "bisyrmsd_angstrom": rmsd,
                "lddt_pli": lddt,
                "ligand_coverage": coverage,
                "rmsd_success": rmsd_success,
                "lddt_pli_success": lddt_success,
                "combined_success": combined_success,
                "confidence_score": confidence.get(
                    "confidence_score", ""
                ),
                "ptm": confidence.get("ptm", ""),
                "iptm": confidence.get("iptm", ""),
                "ligand_iptm": confidence.get(
                    "ligand_iptm", ""
                ),
                "complex_plddt": confidence.get(
                    "complex_plddt", ""
                ),
                "complex_iplddt": confidence.get(
                    "complex_iplddt", ""
                ),
                "complex_ipde": confidence.get(
                    "complex_ipde", ""
                ),
            }

            rows.append(row)
            model_rows.append(row)

        successful = [
            row["pose"]
            for row in model_rows
            if row["combined_success"]
        ]

        best_rmsd_row = min(
            model_rows,
            key=lambda row: row["bisyrmsd_angstrom"],
        )
        best_lddt_row = max(
            model_rows,
            key=lambda row: row["lddt_pli"],
        )

        summary = {
            "pdb_id": pdb_id,
            "model": model,
            "top1_success": model_rows[0]["combined_success"],
            "best_of_5_success": bool(successful),
            "successful_pose_count": len(successful),
            "successful_poses": ";".join(map(str, successful)),
            "top1_rmsd": model_rows[0]["bisyrmsd_angstrom"],
            "top1_lddt_pli": model_rows[0]["lddt_pli"],
            "best_rmsd": best_rmsd_row["bisyrmsd_angstrom"],
            "best_rmsd_pose": best_rmsd_row["pose"],
            "best_lddt_pli": best_lddt_row["lddt_pli"],
            "best_lddt_pose": best_lddt_row["pose"],
            "mean_rmsd": statistics.mean(
                row["bisyrmsd_angstrom"]
                for row in model_rows
            ),
            "mean_lddt_pli": statistics.mean(
                row["lddt_pli"]
                for row in model_rows
            ),
            "affinity_pred_value": "",
            "affinity_probability_binary": "",
        }

        affinity_path = model_dir / f"affinity_{pdb_id}.json"

        if affinity_path.is_file():
            affinity = load_json(affinity_path)
            summary["affinity_pred_value"] = affinity.get(
                "affinity_pred_value", ""
            )
            summary["affinity_probability_binary"] = affinity.get(
                "affinity_probability_binary", ""
            )

        summaries.append(summary)

    pose_output = complex_dir / f"{pdb_id}_pose_results.csv"
    summary_output = complex_dir / f"{pdb_id}_model_summary.csv"

    with pose_output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=rows[0].keys(),
        )
        writer.writeheader()
        writer.writerows(rows)

    with summary_output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=summaries[0].keys(),
        )
        writer.writeheader()
        writer.writerows(summaries)

    print(
        f"{'Model':8s} {'Top-1':>7s} {'Successes':>10s} "
        f"{'Top1 RMSD':>11s} {'Top1 lDDT':>11s} "
        f"{'Best RMSD':>11s} {'Best lDDT':>11s}"
    )

    for summary in summaries:
        print(
            f"{summary['model']:8s} "
            f"{str(summary['top1_success']):>7s} "
            f"{summary['successful_pose_count']:>5d}/5    "
            f"{summary['top1_rmsd']:11.3f} "
            f"{summary['top1_lddt_pli']:11.4f} "
            f"{summary['best_rmsd']:11.3f} "
            f"{summary['best_lddt_pli']:11.4f}"
        )

    print()
    print("Written:", pose_output)
    print("Written:", summary_output)


if __name__ == "__main__":
    main()
