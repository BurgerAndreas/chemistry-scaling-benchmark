#!/usr/bin/env python3
"""Generate small molecule XYZ files (3-20 atoms) for benchmarking."""

from rdkit import Chem
from rdkit.Chem import AllChem
import os

# Define molecules with their SMILES and expected atom count
MOLECULES = [
    # 3 atoms
    (3, "O", "water"),
    (3, "[H]S[H]", "hydrogen_sulfide"),
    (3, "O=C=O", "carbon_dioxide"),
    # 4 atoms
    (4, "N", "ammonia"),
    (4, "P", "phosphine"),
    (4, "[As]", "arsine"),
    # 5 atoms
    (5, "C", "methane"),
    (5, "Si", "silane"),
    (5, "[B-]([H])([H])([H])[H]", "borohydride"),
    # 6 atoms
    (6, "CO", "methanol"),
    (6, "N#N", "dinitrogen"),
    (6, "B2H6", "diborane"),
    # 7 atoms
    (7, "CN", "methylamine"),
    (7, "C1=CC1", "cyclopropene"),
    (7, "C=C=C", "allene"),
    # 8 atoms
    (8, "CC", "ethane"),
    (8, "C1=CC=C1", "cyclobutadiene"),
    (8, "CC#N", "acetonitrile"),
    # 9 atoms
    (9, "CCO", "ethanol"),
    (9, "C1CC1", "cyclopropane"),
    (9, "C=C(C)C", "isobutylene"),
    # 10 atoms
    (10, "CCC", "propane"),
    (10, "C1=CNC=C1", "pyrrole"),
    (10, "C(F)(F)F", "fluoroform"),
    # 11 atoms
    (11, "CC(C)C", "isobutane"),
    (11, "C1=CC=C[CH2]1", "cyclopentadiene"),
    (11, "C1=CC=N1", "pyridine"),
    # 12 atoms
    (12, "C1=CC=CC=C1", "benzene"),
    (12, "C(C)(C)C", "neopentane"),
    (12, "C1CCCCC1", "cyclohexane"),
    # 13 atoms
    (13, "CC1=CC=CC=C1", "toluene"),
    (13, "C1=CC=C(C=C1)O", "phenol"),
    (13, "c1c(N)cccc1", "aniline"),
    # 14 atoms
    (14, "C1=CC=C(C=C1)C=O", "benzaldehyde"),
    (14, "C1=CC2=CC=CC=C2C=C1", "naphthalene"),
    (14, "CC(=O)C1=CC=CC=C1", "acetophenone"),
    # 15 atoms
    (15, "C(C(C(C(C(C=O)O)O)O)O)O", "glucose"),
    (15, "C1=CC=C(C=C1)[C@H](C)O", "phenylethanol"),
    (15, "c1(C(=O)O)ccccc1", "benzoic_acid"),
    # 16 atoms
    (16, "C1=CC=C(C=C1)C(C)(C)C", "tert_butylbenzene"),
    (16, "c1(O)c(O)cccc1", "catechol"),
    (16, "C1=CC=C(C=C1)C#N", "benzonitrile"),
    # 17 atoms
    (17, "C(C(C(=O)O)N)S", "cysteine"),
    (17, "c1(N)c(Cl)cccc1", "2_chloroaniline"),
    (17, "c1(C)c(C)cccc1", "o_xylene"),
    # 18 atoms
    (18, "CC(C)(C)OC(C)(C)C", "di_tert_butyl_ether"),
    (18, "C1=CC=C(C=C1)C(C)=O", "propiophenone"),
    (18, "c1(OC)ccccc1", "anisole"),
    # 19 atoms
    (19, "C(CC(C(=O)O)N)S", "homocysteine"),
    (19, "c1(C)c(C)c(C)ccc1", "1_2_3_trimethylbenzene"),
    (19, "c1ccc(C(c2ccccc2)c2ccccc2)cc1", "triphenylmethane"),
    # 20 atoms
    (20, "C1=CC=C(C=C1)C(=O)C1=CC=CC=C1", "benzophenone"),
    (20, "C1=CC=C(C=C1)C(C)(C)C1=CC=CC=C1", "diphenylmethane"),
    (20, "c1c(Cl)c(Cl)c(Cl)c(Cl)c1Cl", "pentachlorobenzene"),
]


def generate_xyz(smiles, name):
    """Generate XYZ file from SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Failed to create molecule from SMILES: {smiles}")
        return None, None

    mol = Chem.AddHs(mol)
    natoms = mol.GetNumAtoms()

    result = AllChem.EmbedMolecule(mol, randomSeed=42)
    if result != 0:
        print(f"Warning: EmbedMolecule failed for {name}, trying with random coords")
        try:
            AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)
        except:
            print(f"Failed to generate coordinates for {name}")
            return None, None

    try:
        AllChem.UFFOptimizeMolecule(mol)
    except:
        print(f"Warning: UFF optimization failed for {name}")

    conf = mol.GetConformer()
    xyz_lines = [f"{natoms}"]
    xyz_lines.append(f"{name} - generated from SMILES: {smiles}")

    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        symbol = atom.GetSymbol()
        xyz_lines.append(f"{symbol:2s} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}")

    return "\n".join(xyz_lines), natoms


def main():
    output_dir = "xyz_structures"
    os.makedirs(output_dir, exist_ok=True)

    counts = {}

    for _, smiles, name in MOLECULES:
        print(f"Generating {name}...")
        xyz_content, natoms = generate_xyz(smiles, name)

        if xyz_content and natoms:
            # get a unique filename
            count = counts.get(natoms, 0) + 1
            counts[natoms] = count
            filename = f"molecule_{natoms:03d}atoms_{count}.xyz"
            filepath = os.path.join(output_dir, filename)

            with open(filepath, "w") as f:
                f.write(xyz_content)
            print(f"  ✓ Created {filepath}")
        else:
            print(f"  ✗ Failed to create molecule for {name} ({smiles})")


if __name__ == "__main__":
    main()
