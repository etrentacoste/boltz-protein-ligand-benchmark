#!/usr/bin/env python3

import argparse
import csv
import hashlib
import json
import re
import shutil
import urllib.request
from pathlib import Path

import gemmi
from rdkit import Chem


def clean_value(value):
    value = str(value).strip()

    if value in {"", ".", "?"}:
        return ""

    if (
        len(value) >= 2
        and value[0] == value[-1]
        and value[0] in {"'", '"'}
    ):
        value = value[1:-1]

    if value.startswith(";") and value.endswith(";"):
        value = value[1:-1].strip()

    return value.strip()


def normalize_sequence(sequence):
    sequence = "".join(
        clean_value(sequence).split()
    ).upper()

    if not sequence:
        raise RuntimeError("Empty protein sequence")

    if not re.fullmatch(r"[A-Z]+", sequence):
        raise RuntimeError(
            f"Unsupported sequence characters: {sequence}"
        )

    return sequence


def polymer_chain_sequences(block):
    table = block.find(
        "_entity_poly.",
        [
            "entity_id",
            "type",
            "pdbx_seq_one_letter_code_can",
            "pdbx_strand_id",
        ],
    )

    chain_sequences = {}

    for row in table:
        entity_id = clean_value(row[0])
        polymer_type = clean_value(row[1])
        sequence = clean_value(row[2])
        strands = clean_value(row[3])

        if "polypeptide" not in polymer_type.lower():
            continue

        sequence = normalize_sequence(sequence)

        chains = [
            chain.strip()
            for chain in strands.split(",")
            if chain.strip()
        ]

        for chain in chains:
            if (
                chain in chain_sequences
                and chain_sequences[chain]["sequence"]
                != sequence
            ):
                raise RuntimeError(
                    f"Chain {chain} maps to multiple sequences"
                )

            chain_sequences[chain] = {
                "entity_id": entity_id,
                "sequence": sequence,
            }

    return chain_sequences


def sequence_hash(sequence):
    return hashlib.sha256(
        sequence.encode("utf-8")
    ).hexdigest()[:16]


def valid_ccd(path):
    if not path.is_file() or path.stat().st_size < 100:
        return False

    try:
        gemmi.cif.read_file(str(path))
        return True
    except Exception:
        return False


def download_ccd(ccd, output):
    if valid_ccd(output):
        return "existing"

    url = (
        "https://files.rcsb.org/ligands/download/"
        f"{ccd}.cif"
    )
    temporary = output.with_suffix(".cif.part")

    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": (
                "boltz-benchmark-preparation/1.0"
            ),
        },
    )

    with urllib.request.urlopen(
        request,
        timeout=120,
    ) as response:
        temporary.write_bytes(response.read())

    if not valid_ccd(temporary):
        if temporary.exists():
            temporary.unlink()

        raise RuntimeError(
            f"Invalid CCD download for {ccd}"
        )

    temporary.replace(output)
    return "downloaded"


def extract_smiles(ccd_path):
    block = gemmi.cif.read_file(
        str(ccd_path)
    ).sole_block()

    table = block.find(
        "_pdbx_chem_comp_descriptor.",
        [
            "type",
            "program",
            "descriptor",
        ],
    )

    candidates = []

    for row in table:
        descriptor_type = clean_value(row[0])
        program = clean_value(row[1])
        descriptor = clean_value(row[2])

        if "SMILES" not in descriptor_type.upper():
            continue

        if not descriptor:
            continue

        canonical = (
            "CANONICAL" in descriptor_type.upper()
        )

        program_upper = program.upper()

        if canonical and "OPENEYE" in program_upper:
            priority = 0
        elif canonical and "CACTVS" in program_upper:
            priority = 1
        elif canonical:
            priority = 2
        elif "OPENEYE" in program_upper:
            priority = 3
        elif "CACTVS" in program_upper:
            priority = 4
        else:
            priority = 5

        molecule = Chem.MolFromSmiles(descriptor)

        if molecule is None:
            continue

        candidates.append({
            "priority": priority,
            "smiles": descriptor,
            "program": program,
            "type": descriptor_type,
            "heavy_atoms": molecule.GetNumHeavyAtoms(),
            "canonical_rdkit": Chem.MolToSmiles(
                molecule,
                isomericSmiles=True,
            ),
        })

    if not candidates:
        raise RuntimeError(
            f"No RDKit-valid SMILES in {ccd_path}"
        )

    candidates.sort(
        key=lambda item: (
            item["priority"],
            -item["heavy_atoms"],
        )
    )

    return candidates[0], candidates


def split_semicolon(value):
    return [
        item.strip()
        for item in str(value).split(";")
        if item.strip()
    ]


def prepare_one(row, root, ccd_cache):
    pdb_id = row["pdb_id"].upper()
    target_ccd = row["target_ccd"].upper()

    source_cif = (
        root
        / "selection"
        / "full_cifs"
        / f"{pdb_id}.cif"
    )

    if not source_cif.is_file():
        raise FileNotFoundError(source_cif)

    reference_dir = root / "reference" / pdb_id
    reference_dir.mkdir(parents=True, exist_ok=True)

    experimental_cif = (
        reference_dir
        / f"experimental_{pdb_id}.cif"
    )

    shutil.copy2(source_cif, experimental_cif)

    block = gemmi.cif.read_file(
        str(source_cif)
    ).sole_block()

    chain_sequences = polymer_chain_sequences(block)
    contact_chains = split_semicolon(
        row["contact_chains"]
    )

    if not contact_chains:
        raise RuntimeError(
            f"{pdb_id}: no contact chains"
        )

    missing_chains = [
        chain
        for chain in contact_chains
        if chain not in chain_sequences
    ]

    if missing_chains:
        raise RuntimeError(
            f"{pdb_id}: sequence not found for chains "
            f"{missing_chains}; available="
            f"{sorted(chain_sequences)}"
        )

    grouped = {}

    for chain in contact_chains:
        sequence = chain_sequences[chain]["sequence"]
        grouped.setdefault(sequence, []).append(chain)

    proteins = []
    fasta_records = []

    for sequence, chains in grouped.items():
        digest = sequence_hash(sequence)
        msa_path = (
            root
            / "msa"
            / "by_sequence"
            / f"{digest}.csv"
        )

        protein_id = (
            chains[0]
            if len(chains) == 1
            else chains
        )

        proteins.append({
            "id": protein_id,
            "sequence": sequence,
            "msa": str(msa_path),
            "sequence_hash": digest,
            "source_chains": chains,
        })

        fasta_records.append(
            f">{pdb_id}|chains={','.join(chains)}|"
            f"sha256={digest}\n{sequence}\n"
        )

    fasta_path = reference_dir / f"{pdb_id}.fasta"
    fasta_path.write_text(
        "".join(fasta_records)
    )

    ccd_cache.mkdir(parents=True, exist_ok=True)
    cached_ccd = ccd_cache / f"{target_ccd}.cif"
    ccd_status = download_ccd(
        target_ccd,
        cached_ccd,
    )

    local_ccd = reference_dir / f"{target_ccd}.cif"
    shutil.copy2(cached_ccd, local_ccd)

    selected_smiles, all_smiles = extract_smiles(
        local_ccd
    )

    config = {
        "pdb_id": pdb_id,
        "system": row["title"],
        "resolution_angstrom": float(
            row["resolution"]
        ),
        "proteins": proteins,
        "ligands": [
            {
                "id": "L",
                "ccd": target_ccd,
                "smiles": selected_smiles["smiles"],
                "role": "primary_ligand",
                "affinity_binder": True,
            }
        ],
        "reference": {
            "protein_chains": contact_chains,
            "ligand_ccd": target_ccd,
            "ligand_author_chain": (
                row["selected_ligand_chain"]
            ),
            "ligand_residue": (
                row["selected_ligand_residue"]
            ),
            "remove_waters": True,
        },
        "cofactors": [],
        "metals": [],
        "selection_quality": {
            "r_free": float(row["r_free"]),
            "ligand_mean_occupancy": float(
                row["ligand_mean_occupancy"]
            ),
            "ligand_rscc": float(row["ligand_rscc"]),
            "ligand_rsr": float(row["ligand_rsr"]),
            "contact_residue_count": int(
                row["contact_residue_count"]
            ),
            "validation_decision": (
                row["validation_decision"]
            ),
        },
        "smiles_source": {
            "type": selected_smiles["type"],
            "program": selected_smiles["program"],
            "canonical_rdkit": (
                selected_smiles["canonical_rdkit"]
            ),
            "heavy_atoms": (
                selected_smiles["heavy_atoms"]
            ),
        },
        "status": "automatically_prepared",
    }

    config_path = root / "config" / f"{pdb_id}.json"
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(
        json.dumps(config, indent=2) + "\n"
    )

    smiles_audit = (
        reference_dir
        / f"{target_ccd}_smiles_candidates.json"
    )
    smiles_audit.write_text(
        json.dumps(all_smiles, indent=2) + "\n"
    )

    return {
        "pdb_id": pdb_id,
        "target_ccd": target_ccd,
        "contact_chains": ";".join(contact_chains),
        "protein_group_count": len(proteins),
        "sequence_hashes": ";".join(
            protein["sequence_hash"]
            for protein in proteins
        ),
        "selected_smiles": selected_smiles["smiles"],
        "rdkit_canonical_smiles": (
            selected_smiles["canonical_rdkit"]
        ),
        "ligand_heavy_atoms": (
            selected_smiles["heavy_atoms"]
        ),
        "smiles_program": selected_smiles["program"],
        "smiles_type": selected_smiles["type"],
        "ccd_status": ccd_status,
        "status": "SUCCESS",
    }


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--batch",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--root",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--report",
        type=Path,
        required=True,
    )

    args = parser.parse_args()

    with args.batch.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    ccd_cache = args.root / "reference" / "ccd_cache"
    results = []

    for index, row in enumerate(rows, start=1):
        pdb_id = row["pdb_id"].upper()

        try:
            result = prepare_one(
                row,
                args.root,
                ccd_cache,
            )

            print(
                f"[{index:2d}/{len(rows)}] "
                f"{pdb_id}: SUCCESS "
                f"chains={result['contact_chains']} "
                f"sequences={result['sequence_hashes']} "
                f"ligand={result['target_ccd']}"
            )

        except Exception as error:
            result = {
                "pdb_id": pdb_id,
                "status": "FAILED",
                "error": str(error),
            }

            print(
                f"[{index:2d}/{len(rows)}] "
                f"{pdb_id}: FAILED {error}"
            )

        results.append(result)

    fieldnames = []

    for result in results:
        for field in result:
            if field not in fieldnames:
                fieldnames.append(field)

    args.report.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with args.report.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(results)

    failures = [
        result
        for result in results
        if result["status"] != "SUCCESS"
    ]

    unique_sequences = set()

    for result in results:
        if result["status"] == "SUCCESS":
            unique_sequences.update(
                result["sequence_hashes"].split(";")
            )

    print()
    print("Complexes:", len(results))
    print("Successful:", len(results) - len(failures))
    print("Failed:", len(failures))
    print("Unique protein sequences:", len(unique_sequences))
    print("Report:", args.report)

    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
