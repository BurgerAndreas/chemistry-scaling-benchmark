"""
XYZ file parser for molecular structures.
"""

import numpy as np
from typing import Tuple, List


def parse_xyz(xyz_file: str) -> Tuple[List[str], np.ndarray, int]:
    """
    Parse XYZ file and return atomic symbols, coordinates, and number of atoms.

    Args:
        xyz_file: Path to XYZ file

    Returns:
        Tuple of (symbols, coordinates, natoms) where:
        - symbols: List of atomic symbols
        - coordinates: numpy array of shape (natoms, 3) with coordinates in Angstrom
        - natoms: Number of atoms
    """
    with open(xyz_file, "r") as f:
        lines = f.readlines()

    # First line is number of atoms
    natoms = int(lines[0].strip())

    # Second line is comment (skip)

    # Remaining lines are atomic symbols and coordinates
    symbols = []
    coords = []

    for i in range(2, 2 + natoms):
        parts = lines[i].strip().split()
        symbols.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

    coords = np.array(coords)

    return symbols, coords, natoms


def xyz_to_pyscf_format(symbols: List[str], coords: np.ndarray) -> str:
    """
    Convert symbols and coordinates to PySCF format string.

    Args:
        symbols: List of atomic symbols
        coords: Coordinates array (natoms, 3) in Angstrom

    Returns:
        PySCF molecule string in format "SYMBOL X Y Z; SYMBOL X Y Z; ..."
    """
    mol_str_parts = []
    for i, symbol in enumerate(symbols):
        mol_str_parts.append(f"{symbol} {coords[i, 0]} {coords[i, 1]} {coords[i, 2]}")

    return "; ".join(mol_str_parts)
