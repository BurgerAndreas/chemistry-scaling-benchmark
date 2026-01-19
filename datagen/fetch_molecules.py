#!/usr/bin/env python3
"""Fetch XYZ geometries of molecules with specific atom counts from PubChem."""

import os
import requests
import time
from rdkit import Chem
from rdkit.Chem import AllChem

# Target atom counts
TARGET_ATOMS = list(range(10, 101, 10))

# Create output directory
os.makedirs("xyz_structures", exist_ok=True)


def get_xyz_from_mol(mol, name="molecule"):
    """Convert RDKit mol to XYZ format string."""
    conf = mol.GetConformer()
    num_atoms = mol.GetNumAtoms()

    lines = [str(num_atoms), name]
    for i, atom in enumerate(mol.GetAtoms()):
        pos = conf.GetAtomPosition(i)
        symbol = atom.GetSymbol()
        lines.append(f"{symbol:2s} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}")

    return "\n".join(lines)


def fetch_sdf_by_cid(cid, prefer_3d=True):
    """Fetch SDF structure from PubChem by CID."""
    if prefer_3d:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                return response.text, True
        except:
            pass

    # Fall back to 2D
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF"
    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            return response.text, False
    except:
        pass
    return None, False


def get_molecule_name(cid):
    """Get molecule name from PubChem."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName,Title/JSON"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            props = data["PropertyTable"]["Properties"][0]
            return props.get("Title", props.get("IUPACName", f"CID_{cid}"))
    except:
        pass
    return f"CID_{cid}"


def process_molecule(cid, target_atoms, tolerance=0):
    """Process a molecule and check if it has the target atom count."""
    sdf, is_3d = fetch_sdf_by_cid(cid)
    if sdf is None:
        return None, None, None

    mol = Chem.MolFromMolBlock(sdf)
    if mol is None:
        return None, None, None

    # Add hydrogens
    mol = Chem.AddHs(mol)

    # Generate 3D coordinates if needed
    try:
        if not is_3d or mol.GetNumConformers() == 0:
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except:
        try:
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        except:
            return None, None, None

    if mol.GetNumConformers() == 0:
        return None, None, None

    num_atoms = mol.GetNumAtoms()

    if abs(num_atoms - target_atoms) <= tolerance:
        return mol, num_atoms, get_molecule_name(cid)

    return None, num_atoms, None


def get_cids_by_heavy_atoms(heavy_count, limit=200):
    """Get CIDs for compounds with specific heavy atom count using property table."""
    # Use the listkey approach for large searches
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/1-{50000 + heavy_count * 2000}/property/HeavyAtomCount/JSON"

    # Alternative: sample from known CID ranges based on molecule size
    # Small molecules (low CIDs), larger molecules (higher CIDs)
    return []


# Curated molecules by approximate total atom count (including H)
# Format: (CID, common_name, approx_total_atoms)
CANDIDATE_MOLECULES = {
    10: [
        (6342, "propanol", 12),  # C3H8O
        (702, "ethanol", 9),  # C2H6O
        (887, "methanol", 6),  # CH4O
        (6344, "butanol", 15),  # C4H10O
        (180, "acetone", 10),  # C3H6O
        (679, "dimethyl ether", 9),  # C2H6O
        (6569, "propane", 11),  # C3H8
        (8857, "cyclopropane", 9),  # C3H6
        (6334, "propane", 11),  # C3H8
        (7843, "propene", 9),  # C3H6
    ],
    20: [
        (241, "benzene", 12),  # C6H6
        (7500, "cyclohexane", 18),  # C6H12
        (6212, "hexane", 20),  # C6H14
        (7843, "propene", 9),
        (8058, "styrene", 16),  # C8H8
        (7501, "cyclohexene", 16),  # C6H10
        (996, "phenol", 13),  # C6H6O
        (243, "benzoic acid", 15),  # C7H6O2
        (1140, "aniline", 14),  # C6H7N
        (7809, "toluene", 15),  # C7H8
        (31276, "heptane", 23),  # C7H16
        (356, "octane", 26),  # C8H18
    ],
    30: [
        (931, "naphthalene", 18),  # C10H8
        (6549, "decane", 32),  # C10H22
        (2519, "caffeine", 24),  # C8H10N4O2
        (2244, "aspirin", 21),  # C9H8O4
        (1983, "acetaminophen", 20),  # C8H9NO2
        (3386, "ibuprofen", 33),  # C13H18O2
        (1032, "pyridine", 11),  # C5H5N
        (2381, "dopamine", 22),  # C8H11NO2
        (3348, "quercetin", 32),  # C15H10O7
        (1423, "nicotinic acid", 14),  # C6H5NO2
        (5950, "adrenaline", 23),  # C9H13NO3
    ],
    40: [
        (2519, "caffeine", 24),
        (3386, "ibuprofen", 33),
        (2662, "kaempferol", 29),  # C15H10O6
        (5245, "naproxen", 31),  # C14H14O3
        (2157, "glucose", 24),  # C6H12O6
        (5988, "sucrose", 45),  # C12H22O11
        (5754, "cortisone", 44),  # C21H28O5
        (3033, "lidocaine", 46),  # C14H22N2O
        (71587, "donepezil", 50),  # C24H29NO3
        (5743, "morphine", 40),  # C17H19NO3
    ],
    50: [
        (5743, "morphine", 40),
        (5288, "codeine", 43),  # C18H21NO3
        (5284371, "paclitaxel", 113),
        (5904, "penicillin G", 41),  # C16H18N2O4S
        (5994, "progesterone", 48),  # C21H30O2
        (5757, "estradiol", 44),  # C18H24O2
        (6013, "testosterone", 49),  # C19H28O2
        (11979, "melatonin", 33),  # C13H16N2O2
        (3034034, "tamoxifen", 54),  # C26H29NO
        (2244, "aspirin", 21),
        (5281515, "apigenin", 28),  # C15H10O5
        (4485, "methotrexate", 54),  # C20H22N8O5
    ],
    60: [
        (5757, "estradiol", 44),
        (5994, "progesterone", 48),
        (5997, "cholesterol", 74),  # C27H46O
        (3121, "fluoxetine", 43),  # C17H18F3NO
        (2554, "chloroquine", 48),  # C18H26ClN3
        (60823, "simvastatin", 67),  # C25H38O5
        (54454, "losartan", 57),  # C22H23ClN6O
        (3676, "metformin", 17),  # C4H11N5
        (2082, "clonidine", 21),  # C9H9Cl2N3
        (54677470, "rivaroxaban", 57),  # C19H18ClN3O5S
        (60606, "sildenafil", 63),  # C22H30N6O4S
    ],
    70: [
        (5997, "cholesterol", 74),
        (2244, "atorvastatin", 83),  # C33H35FN2O5 (CID 60823 is wrong)
        (60823, "simvastatin", 67),
        (154234, "ivermectin", 157),  # C48H74O14
        (16078, "thyroxine", 27),  # C15H11I4NO4
        (5280343, "quercetin", 32),
        (6918837, "sunitinib", 58),  # C22H27FN4O2
        (2435, "diclofenac", 31),  # C14H11Cl2NO2
        (54677470, "rivaroxaban", 57),
        (3062316, "gefitinib", 56),  # C22H24ClFN4O3
        (71617, "nifedipine", 51),  # C17H18N2O6
    ],
    80: [
        (5997, "cholesterol", 74),
        (446220, "ubiquinone", 133),  # C59H90O4
        (5280961, "genistein", 28),  # C15H10O5
        (6251, "isoleucine", 22),  # C6H13NO2
        (6106, "leucine", 22),  # C6H13NO2
        (5951, "serotonin", 23),  # C10H12N2O
        (3821, "ketoconazole", 63),  # C26H28Cl2N4O4
        (68617, "nebivolol", 74),  # C22H25F2NO4
        (3002977, "fullerene C60", 60),  # C60
    ],
    90: [
        (5997, "cholesterol", 74),
        (446220, "ubiquinone", 133),
        (11110, "biotin", 24),  # C10H16N2O3S
        (5280793, "tocopherol", 81),  # C29H50O2
        (3002977, "fullerene C60", 60),
        (11954189, "vitamin D3", 72),  # C27H44O
        (5281, "retinol", 49),  # C20H30O
        (5280794, "tocotrienol", 79),  # C29H44O2
    ],
    100: [
        (5280793, "tocopherol", 81),
        (446220, "ubiquinone", 133),
        (5280794, "tocotrienol", 79),
        (6857552, "calcifediol", 66),  # C27H44O2
        (73581, "vancomycin aglycon", 129),  # large
        (5282411, "astaxanthin", 82),  # C40H52O4
        (5281224, "lycopene", 96),  # C40H56
        (5280489, "beta-carotene", 96),  # C40H56
    ],
}


def find_exact_match(target, tolerance=0):
    """Search through candidate molecules to find one with exact atom count."""
    print(f"  Searching candidates for {target} atoms (tolerance ±{tolerance})...")

    # Also check nearby targets for potential matches
    check_targets = [target]
    if tolerance > 0:
        check_targets.extend([target - 10, target + 10])

    for check_target in check_targets:
        candidates = CANDIDATE_MOLECULES.get(check_target, [])
        for cid, name, approx in candidates:
            print(
                f"    Trying {name} (CID {cid}, ~{approx} atoms)...",
                end=" ",
                flush=True,
            )

            mol, num_atoms, mol_name = process_molecule(cid, target, tolerance)

            if mol is not None:
                print(f"MATCH! {num_atoms} atoms")
                return mol, cid, mol_name, num_atoms
            elif num_atoms is not None:
                print(f"{num_atoms} atoms (not matching)")
            else:
                print("failed to process")

            time.sleep(0.2)

    return None, None, None, None


def scan_cid_range(target, start_cid, end_cid, step=50, tolerance=0):
    """Scan a CID range to find molecules with target atom count."""
    print(f"  Scanning CID range {start_cid}-{end_cid}...")

    for cid in range(start_cid, end_cid, step):
        mol, num_atoms, name = process_molecule(cid, target, tolerance)

        if mol is not None:
            return mol, cid, name, num_atoms

        time.sleep(0.15)

    return None, None, None, None


def main():
    results = {}

    for target in TARGET_ATOMS:
        print(f"\n{'=' * 60}")
        print(f"TARGET: {target} atoms")
        print("=" * 60)

        # First try curated candidates with exact match
        mol, cid, name, num_atoms = find_exact_match(target, tolerance=0)

        # If no exact match, try with tolerance
        if mol is None:
            print(f"\n  No exact match. Trying with ±1 tolerance...")
            mol, cid, name, num_atoms = find_exact_match(target, tolerance=1)

        # If still no match, scan CID ranges based on target size
        if mol is None:
            print(f"\n  Scanning additional CID ranges...")
            if target <= 30:
                ranges = [(100, 2000, 20), (2000, 10000, 50)]
            elif target <= 60:
                ranges = [(5000, 30000, 100), (30000, 80000, 200)]
            else:
                ranges = [(50000, 200000, 300), (200000, 500000, 500)]

            for start, end, step in ranges:
                mol, cid, name, num_atoms = scan_cid_range(
                    target, start, end, step, tolerance=1
                )
                if mol is not None:
                    break

        if mol is not None:
            xyz_content = get_xyz_from_mol(mol, f"{name} (CID {cid})")
            filename = f"xyz_structures/molecule_{target:03d}atoms.xyz"
            with open(filename, "w") as f:
                f.write(xyz_content)
            print(f"\n  SAVED: {filename} ({num_atoms} atoms)")
            results[target] = (name, cid, num_atoms, filename)
        else:
            print(f"\n  FAILED: Could not find molecule with ~{target} atoms")
            results[target] = None

        time.sleep(0.5)

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY OF SAVED XYZ FILES")
    print("=" * 70)
    print(f"{'Target':>6} | {'Actual':>6} | {'Name':<40} | CID")
    print("-" * 70)
    for target in TARGET_ATOMS:
        if results.get(target):
            name, cid, num_atoms, filename = results[target]
            print(f"{target:>6} | {num_atoms:>6} | {name[:40]:<40} | {cid}")
        else:
            print(f"{target:>6} | {'--':>6} | {'NOT FOUND':<40} | --")

    print("\nFiles saved in: xyz_structures/")


if __name__ == "__main__":
    main()
