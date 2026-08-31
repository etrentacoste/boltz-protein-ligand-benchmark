import argparse
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D


BOND_TYPES = {
    1: Chem.BondType.SINGLE,
    2: Chem.BondType.DOUBLE,
    3: Chem.BondType.TRIPLE,
    4: Chem.BondType.AROMATIC,
}


def decode_text(value):
    """Decode strings used by different Boltz NPZ formats."""
    if isinstance(value, str):
        return value.strip()

    if isinstance(value, bytes):
        return value.decode().strip()

    if np.isscalar(value):
        return str(value).strip()

    characters = []

    for number in value:
        number = int(number)
        if number != 0:
            characters.append(chr(number + 32))

    return "".join(characters).strip()


def read_ligand_coordinates(cif_path, ligand_chain):
    structure = gemmi.read_structure(str(cif_path))
    coordinates = {}

    for model in structure:
        for chain in model:
            for residue in chain:
                is_ligand = (
                    chain.name == ligand_chain
                    or residue.name.startswith("LIG")
                )

                if not is_ligand:
                    continue

                for atom in residue:
                    if atom.element.name == "H":
                        continue

                    name = atom.name.strip()

                    if name in coordinates:
                        raise RuntimeError(
                            f"Duplicate ligand atom name {name} "
                            f"in {cif_path}"
                        )

                    coordinates[name] = (
                        atom.pos.x,
                        atom.pos.y,
                        atom.pos.z,
                        atom.element.name,
                    )

    return coordinates


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Create explicit ligand SDF files from Boltz NPZ topology "
            "and predicted mmCIF coordinates."
        )
    )

    parser.add_argument("--pdb-id", required=True)
    parser.add_argument(
        "--model",
        required=True,
        choices=["boltz1", "boltz2"],
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path.home() / "boltz_visualization",
    )
    parser.add_argument("--ligand-chain", default="L")
    parser.add_argument("--poses", type=int, default=5)

    args = parser.parse_args()

    pdb_id = args.pdb_id.upper()
    model_dir = args.root / pdb_id / "production" / args.model
    npz_path = model_dir / f"{pdb_id}_processed.npz"
    sdf_dir = model_dir / "ligand_sdf"

    if not npz_path.is_file():
        raise FileNotFoundError(npz_path)

    data = np.load(npz_path, allow_pickle=True)

    atoms = data["atoms"]
    bonds = data["bonds"]
    chains = data["chains"]

    ligand_chain_record = None

    for chain in chains:
        chain_name = decode_text(chain["name"])

        if chain_name == args.ligand_chain:
            ligand_chain_record = chain
            break

    if ligand_chain_record is None:
        available = [decode_text(chain["name"]) for chain in chains]
        raise RuntimeError(
            f"Ligand chain {args.ligand_chain!r} was not found. "
            f"Available chains: {available}"
        )

    start = int(ligand_chain_record["atom_idx"])
    count = int(ligand_chain_record["atom_num"])
    end = start + count

    ligand_names = [
        decode_text(value)
        for value in atoms["name"][start:end]
    ]

    if len(set(ligand_names)) != len(ligand_names):
        raise RuntimeError(
            f"Ligand atom names are not unique: {ligand_names}"
        )

    ligand_bonds = [
        bond
        for bond in bonds
        if (
            start <= int(bond["atom_1"]) < end
            and start <= int(bond["atom_2"]) < end
        )
    ]

    print("PDB:", pdb_id)
    print("Model:", args.model)
    print("Ligand chain:", args.ligand_chain)
    print("Ligand atom range:", start, "to", end - 1)
    print("Ligand atoms:", len(ligand_names))
    print("Ligand bonds:", len(ligand_bonds))
    print("Atom names:", ligand_names)

    sdf_dir.mkdir(parents=True, exist_ok=True)

    for pose_index in range(args.poses):
        cif_path = model_dir / f"{pdb_id}_model_{pose_index}.cif"
        output_path = (
            sdf_dir
            / f"{pdb_id}_model_{pose_index}_ligand.sdf"
        )

        if not cif_path.is_file():
            raise FileNotFoundError(cif_path)

        coordinates = read_ligand_coordinates(
            cif_path,
            args.ligand_chain,
        )

        missing = [
            name for name in ligand_names
            if name not in coordinates
        ]

        if missing:
            raise RuntimeError(
                f"Atoms missing from {cif_path.name}: {missing}"
            )

        editable = Chem.RWMol()

        for local_index, name in enumerate(ligand_names):
            element = coordinates[name][3]
            atom = Chem.Atom(element)
            atom.SetProp("_TriposAtomName", name)

            if "charge" in atoms.dtype.names:
                charge = int(atoms["charge"][start + local_index])
                atom.SetFormalCharge(charge)

            editable.AddAtom(atom)

        for bond in ligand_bonds:
            atom_1 = int(bond["atom_1"]) - start
            atom_2 = int(bond["atom_2"]) - start
            numeric_type = int(bond["type"])

            if numeric_type not in BOND_TYPES:
                raise RuntimeError(
                    f"Unknown bond type: {numeric_type}"
                )

            editable.AddBond(
                atom_1,
                atom_2,
                BOND_TYPES[numeric_type],
            )

            if numeric_type == 4:
                editable.GetAtomWithIdx(atom_1).SetIsAromatic(True)
                editable.GetAtomWithIdx(atom_2).SetIsAromatic(True)

        molecule = editable.GetMol()
        molecule.SetProp(
            "_Name",
            f"{pdb_id}_{args.model}_model_{pose_index}",
        )

        conformer = Chem.Conformer(count)
        conformer.Set3D(True)

        for atom_index, name in enumerate(ligand_names):
            x, y, z, _ = coordinates[name]
            conformer.SetAtomPosition(
                atom_index,
                Point3D(x, y, z),
            )

        molecule.AddConformer(conformer, assignId=True)

        try:
            Chem.SanitizeMol(molecule)
        except Exception as error:
            print(
                f"Warning: RDKit sanitization for "
                f"pose {pose_index}: {error}"
            )

        writer = Chem.SDWriter(str(output_path))
        writer.SetKekulize(False)
        # OpenStructure's SDF reader does not accept MDL bond type 4.
        # Preserve atom connectivity and coordinates, but encode aromatic
        # bonds as ordinary single bonds in this auxiliary topology file.
        for bond in molecule.GetBonds():
            if (
                bond.GetIsAromatic()
                or bond.GetBondType() == Chem.BondType.AROMATIC
            ):
                bond.SetIsAromatic(False)
                bond.SetBondType(Chem.BondType.SINGLE)

        for atom in molecule.GetAtoms():
            atom.SetIsAromatic(False)

        writer.write(molecule)
        writer.close()

        check = Chem.SDMolSupplier(
            str(output_path),
            removeHs=False,
            sanitize=False,
        )[0]

        if check is None:
            raise RuntimeError(
                f"Could not read generated SDF: {output_path}"
            )

        if check.GetNumAtoms() != count:
            raise RuntimeError(
                f"Incorrect atom count in {output_path}"
            )

        if check.GetNumBonds() != len(ligand_bonds):
            raise RuntimeError(
                f"Incorrect bond count in {output_path}"
            )

        print(
            output_path.name,
            "atoms =", check.GetNumAtoms(),
            "bonds =", check.GetNumBonds(),
        )

    print("SDF generation: SUCCESS")


if __name__ == "__main__":
    main()
