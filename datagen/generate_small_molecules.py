#!/usr/bin/env python3
"""Generate small molecule XYZ files (3-9 atoms) for benchmarking."""

from rdkit import Chem
from rdkit.Chem import AllChem
import os

# Define molecules with their SMILES and expected atom count
MOLECULES = [
    (3, "O", "water"),  # H2O
    (4, "N", "ammonia"),  # NH3
    (5, "C", "methane"),  # CH4
    (6, "CO", "methanol"),  # CH3OH
    (7, "CN", "methylamine"),  # CH3NH2
    (8, "CC", "ethane"),  # C2H6
    (9, "CCO", "ethanol"),  # C2H5OH
]


def generate_xyz(smiles, natoms, name):
    """Generate XYZ file from SMILES string."""
    # Create molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Failed to create molecule from SMILES: {smiles}")
        return None

    # Add hydrogens
    mol = Chem.AddHs(mol)

    # Verify atom count
    actual_natoms = mol.GetNumAtoms()
    if actual_natoms != natoms:
        print(f"Warning: {name} has {actual_natoms} atoms, expected {natoms}")
        natoms = actual_natoms

    # Generate 3D coordinates
    result = AllChem.EmbedMolecule(mol, randomSeed=42)
    if result != 0:
        print(f"Warning: EmbedMolecule failed for {name}, trying with random coords")
        AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)

    # Optimize geometry with UFF
    try:
        AllChem.UFFOptimizeMolecule(mol)
    except:
        print(f"Warning: UFF optimization failed for {name}")

    # Create XYZ format
    conf = mol.GetConformer()
    xyz_lines = [f"{natoms}"]
    xyz_lines.append(f"{name} - generated from SMILES: {smiles}")

    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        symbol = atom.GetSymbol()
        xyz_lines.append(f"{symbol:2s} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}")

    return "\n".join(xyz_lines)


def main():
    output_dir = "xyz_structures"
    os.makedirs(output_dir, exist_ok=True)

    for natoms, smiles, name in MOLECULES:
        filename = f"molecule_{natoms:03d}atoms.xyz"
        filepath = os.path.join(output_dir, filename)

        if os.path.exists(filepath):
            print(f"Skipping {filename} (already exists)")
            continue

        print(f"Generating {filename} ({name})...")
        xyz_content = generate_xyz(smiles, natoms, name)

        if xyz_content:
            with open(filepath, "w") as f:
                f.write(xyz_content)
            print(f"  ✓ Created {filename}")
        else:
            print(f"  ✗ Failed to create {filename}")


if __name__ == "__main__":
    main()
