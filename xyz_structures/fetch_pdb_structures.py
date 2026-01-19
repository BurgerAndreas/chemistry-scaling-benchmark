#!/usr/bin/env python3
"""
Fetch XYZ geometries from RCSB Protein Data Bank for structures with 100-1000 atoms.

Strategy: Search by polymer residue count (sequence length) since:
- deposited_atom_count includes waters and heteroatoms
- Average ~7.5 heavy atoms per amino acid residue
- With hydrogens, ~15 atoms per residue total

For target of N atoms, search for ~N/15 to ~N/7 residues.
"""

import os
import requests
import time
import gzip
from io import StringIO

# Target atom counts from 100 to 1000 in steps of 10
TARGET_ATOMS = list(range(100, 1001, 10))

# Output directory
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))


def search_pdb_by_residue_count(min_res, max_res, limit=20):
    """
    Search RCSB PDB for structures with residue count in given range.
    """
    url = "https://search.rcsb.org/rcsbsearch/v2/query"

    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.deposited_polymer_monomer_count",
                        "operator": "range",
                        "value": {
                            "from": min_res,
                            "to": max_res,
                            "include_lower": True,
                            "include_upper": True,
                        },
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.polymer_entity_count_protein",
                        "operator": "equals",
                        "value": 1,
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": limit},
            "sort": [
                {"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}
            ],
        },
    }

    try:
        response = requests.post(url, json=query, timeout=30)
        if response.status_code == 200:
            data = response.json()
            return [r["identifier"] for r in data.get("result_set", [])]
    except Exception as e:
        print(f"    Search error: {e}")

    return []


def get_pdb_info(pdb_id):
    """Get information about a PDB entry."""
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    try:
        response = requests.get(url, timeout=15)
        if response.status_code == 200:
            return response.json()
    except:
        pass
    return None


def download_pdb_file(pdb_id):
    """Download PDB file from RCSB."""
    # Try regular PDB first
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        response = requests.get(url, timeout=60)
        if response.status_code == 200:
            return response.text
    except:
        pass

    # Try compressed
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb.gz"
    try:
        response = requests.get(url, timeout=60)
        if response.status_code == 200:
            return gzip.decompress(response.content).decode("utf-8")
    except:
        pass

    return None


def pdb_to_xyz(pdb_content, title="molecule", include_hydrogens=True):
    """
    Convert PDB format to XYZ format.
    Only extracts ATOM records (polymer atoms, not waters/ligands).
    """
    atoms = []

    for line in pdb_content.split("\n"):
        # Only use ATOM records (polymer), skip HETATM (waters, ligands)
        if not line.startswith("ATOM"):
            continue

        try:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())

            # Get element symbol
            element = line[76:78].strip()
            if not element:
                atom_name = line[12:16].strip()
                element = "".join(c for c in atom_name[:2] if c.isalpha())
                if len(element) > 1:
                    element = element[0].upper() + element[1].lower()
                else:
                    element = element.upper()

            # Skip hydrogens if not wanted
            if not include_hydrogens and element == "H":
                continue

            if element:
                atoms.append((element, x, y, z))
        except (ValueError, IndexError):
            continue

    if not atoms:
        return None

    lines = [str(len(atoms)), title]
    for elem, x, y, z in atoms:
        lines.append(f"{elem:2s} {x:12.6f} {y:12.6f} {z:12.6f}")

    return "\n".join(lines)


def find_structure_for_target(target, tolerance=10):
    """Find a PDB structure with approximately the target number of atoms."""
    print(f"\nSearching for structure with ~{target} atoms...")

    # Estimate residue count: ~7-8 heavy atoms per residue
    # We search by residues since that's more reliable
    # Most PDB structures don't have hydrogens, so use ~7 atoms/residue
    est_residues = target // 7

    # Search with different residue ranges
    for res_offset in [0, 2, 5, 10, 15, 20]:
        min_res = max(1, est_residues - res_offset)
        max_res = est_residues + res_offset

        print(f"  Residues: {min_res}-{max_res}...", end=" ", flush=True)
        pdb_ids = search_pdb_by_residue_count(min_res, max_res, limit=10)

        if not pdb_ids:
            print("no results")
            time.sleep(0.3)
            continue

        print(f"found {len(pdb_ids)} entries")

        for pdb_id in pdb_ids:
            print(f"    Trying {pdb_id}...", end=" ", flush=True)

            # Get PDB info for title
            info = get_pdb_info(pdb_id)
            if info:
                title = info.get("struct", {}).get("title", pdb_id)[:80]
            else:
                title = pdb_id

            # Download PDB file
            pdb_content = download_pdb_file(pdb_id)
            if not pdb_content:
                print("download failed")
                time.sleep(0.3)
                continue

            # Convert to XYZ (polymer atoms only)
            xyz_content = pdb_to_xyz(pdb_content, f"{title} (PDB: {pdb_id})")
            if not xyz_content:
                print("conversion failed")
                time.sleep(0.3)
                continue

            actual_atoms = int(xyz_content.split("\n")[0])
            print(f"{actual_atoms} atoms")

            if abs(actual_atoms - target) <= tolerance:
                return xyz_content, pdb_id, title, actual_atoms

            time.sleep(0.2)

        time.sleep(0.3)

    return None, None, None, None


def main():
    results = {}
    found_atoms = set()  # Track what we've already found to avoid duplicates

    print("=" * 70)
    print("FETCHING STRUCTURES FROM RCSB PROTEIN DATA BANK")
    print("Target: 100-1000 atoms in steps of 10")
    print("=" * 70)

    for target in TARGET_ATOMS:
        print(f"\n{'=' * 60}")
        print(f"TARGET: {target} atoms")
        print("=" * 60)

        xyz_content, pdb_id, title, actual_atoms = find_structure_for_target(target)

        if xyz_content and actual_atoms not in found_atoms:
            filename = f"molecule_{target:04d}atoms.xyz"
            filepath = os.path.join(OUTPUT_DIR, filename)
            with open(filepath, "w") as f:
                f.write(xyz_content)
            print(f"\n  SAVED: {filename} ({actual_atoms} atoms)")
            results[target] = (pdb_id, title[:50], actual_atoms)
            found_atoms.add(actual_atoms)
        elif xyz_content:
            print(f"\n  SKIPPED: {actual_atoms} atoms already saved")
            results[target] = None
        else:
            print(f"\n  FAILED: Could not find structure with ~{target} atoms")
            results[target] = None

        time.sleep(0.3)

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY OF SAVED XYZ FILES")
    print("=" * 80)
    print(f"{'Target':>6} | {'Actual':>6} | {'PDB ID':<8} | Description")
    print("-" * 80)

    found = 0
    for target in TARGET_ATOMS:
        if results.get(target):
            pdb_id, title, actual_atoms = results[target]
            print(f"{target:>6} | {actual_atoms:>6} | {pdb_id:<8} | {title}")
            found += 1
        else:
            print(f"{target:>6} | {'--':>6} | {'--':<8} | NOT FOUND")

    print(f"\nFound {found}/{len(TARGET_ATOMS)} structures")
    print(f"Files saved in: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
