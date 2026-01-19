#!/usr/bin/env python3
"""Fill gaps in molecule collection by using higher tolerance and different search strategies."""

import os
import requests
import time
import gzip

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

# Find missing targets
ALL_TARGETS = list(range(100, 1001, 10))
existing_files = os.listdir(OUTPUT_DIR)
found_targets = set()
for f in existing_files:
    if f.startswith("molecule_") and f.endswith("atoms.xyz"):
        try:
            target = int(f.split("_")[1].replace("atoms.xyz", ""))
            if 100 <= target <= 1000:
                found_targets.add(target)
        except:
            pass

MISSING_TARGETS = sorted([t for t in ALL_TARGETS if t not in found_targets])
print(f"Missing targets: {len(MISSING_TARGETS)}")
print(MISSING_TARGETS)


def search_pdb_single_chain(min_res, max_res, limit=30):
    """Search for single-chain proteins."""
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
                        "value": {"from": min_res, "to": max_res}
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.polymer_entity_count",
                        "operator": "equals",
                        "value": 1
                    }
                }
            ]
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": limit}}
    }
    try:
        response = requests.post(url, json=query, timeout=30)
        if response.status_code == 200:
            data = response.json()
            return [r["identifier"] for r in data.get("result_set", [])]
    except Exception as e:
        print(f"  Error: {e}")
    return []


def download_pdb(pdb_id):
    """Download PDB file."""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        response = requests.get(url, timeout=60)
        if response.status_code == 200:
            return response.text
    except:
        pass
    return None


def get_pdb_title(pdb_id):
    """Get PDB entry title."""
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    try:
        response = requests.get(url, timeout=15)
        if response.status_code == 200:
            return response.json().get("struct", {}).get("title", pdb_id)[:80]
    except:
        pass
    return pdb_id


def pdb_to_xyz(pdb_content, title="molecule"):
    """Convert PDB to XYZ (only ATOM records)."""
    atoms = []
    for line in pdb_content.split('\n'):
        if not line.startswith('ATOM'):
            continue
        try:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            element = line[76:78].strip()
            if not element:
                atom_name = line[12:16].strip()
                element = ''.join(c for c in atom_name[:2] if c.isalpha())
                if len(element) > 1:
                    element = element[0].upper() + element[1].lower()
                else:
                    element = element.upper()
            if element:
                atoms.append((element, x, y, z))
        except:
            continue

    if not atoms:
        return None

    lines = [str(len(atoms)), title]
    for elem, x, y, z in atoms:
        lines.append(f"{elem:2s} {x:12.6f} {y:12.6f} {z:12.6f}")
    return '\n'.join(lines)


def find_for_target(target, max_tolerance=20):
    """Find structure for target with increasing tolerance."""
    print(f"\nTarget: {target} atoms")

    # Try different tolerances
    for tolerance in [5, 10, 15, 20, 25, 30]:
        if tolerance > max_tolerance:
            break

        est_residues = target // 7

        for res_offset in [0, 3, 6, 10, 15]:
            min_res = max(1, est_residues - res_offset)
            max_res = est_residues + res_offset

            pdb_ids = search_pdb_single_chain(min_res, max_res, limit=20)

            for pdb_id in pdb_ids:
                pdb_content = download_pdb(pdb_id)
                if not pdb_content:
                    continue

                title = get_pdb_title(pdb_id)
                xyz_content = pdb_to_xyz(pdb_content, f"{title} (PDB: {pdb_id})")
                if not xyz_content:
                    continue

                actual_atoms = int(xyz_content.split('\n')[0])

                if abs(actual_atoms - target) <= tolerance:
                    print(f"  Found: {pdb_id} with {actual_atoms} atoms (tolerance: {tolerance})")
                    return xyz_content, pdb_id, title, actual_atoms

                time.sleep(0.1)

            time.sleep(0.2)

    return None, None, None, None


def main():
    results = {}

    for target in MISSING_TARGETS:
        xyz_content, pdb_id, title, actual_atoms = find_for_target(target, max_tolerance=25)

        if xyz_content:
            filename = f"molecule_{target:04d}atoms.xyz"
            filepath = os.path.join(OUTPUT_DIR, filename)
            with open(filepath, "w") as f:
                f.write(xyz_content)
            print(f"  SAVED: {filename} ({actual_atoms} atoms)")
            results[target] = (pdb_id, title[:50], actual_atoms)
        else:
            print(f"  FAILED")
            results[target] = None

        time.sleep(0.3)

    # Summary
    print("\n" + "="*70)
    print("GAP FILLING SUMMARY")
    print("="*70)
    found = sum(1 for v in results.values() if v is not None)
    print(f"Found {found}/{len(MISSING_TARGETS)} missing structures")


if __name__ == "__main__":
    main()
