"""
xTB calculator for GFN2-xTB tight binding method.
"""

from typing import Dict, Any
import numpy as np
from xtb.interface import Calculator
from xtb.utils import get_method
from calculators.base_calculator import BaseCalculator


class XTBCalculator(BaseCalculator):
    """
    Calculator for GFN2-xTB (Tight Binding) using xtb-python.
    """

    def __init__(self):
        super().__init__(
            method_category="tight_binding", method_name="GFN2-xTB", basis_set="N/A"
        )

    def calculate(self, xyz_file: str) -> Dict[str, Any]:
        """
        Perform GFN2-xTB calculation.

        Args:
            xyz_file: Path to XYZ file

        Returns:
            Dictionary with success, energy_hartree, nbasis, error_message
        """
        try:
            # Parse XYZ file
            symbols, coords, natoms = self.parse_xyz(xyz_file)

            # Convert symbols to atomic numbers
            atomic_numbers = []
            atomic_symbol_to_number = {
                "H": 1,
                "He": 2,
                "Li": 3,
                "Be": 4,
                "B": 5,
                "C": 6,
                "N": 7,
                "O": 8,
                "F": 9,
                "Ne": 10,
                "Na": 11,
                "Mg": 12,
                "Al": 13,
                "Si": 14,
                "P": 15,
                "S": 16,
                "Cl": 17,
                "Ar": 18,
                "K": 19,
                "Ca": 20,
                "Sc": 21,
                "Ti": 22,
                "V": 23,
                "Cr": 24,
                "Mn": 25,
                "Fe": 26,
                "Co": 27,
                "Ni": 28,
                "Cu": 29,
                "Zn": 30,
                "Ga": 31,
                "Ge": 32,
                "As": 33,
                "Se": 34,
                "Br": 35,
                "Kr": 36,
            }

            for symbol in symbols:
                if symbol not in atomic_symbol_to_number:
                    return {
                        "success": False,
                        "error_message": f"Unsupported element: {symbol}",
                    }
                atomic_numbers.append(atomic_symbol_to_number[symbol])

            atomic_numbers = np.array(atomic_numbers)

            # Convert coordinates from Angstrom to Bohr (1 Angstrom = 1.88973 Bohr)
            coords_bohr = coords * 1.88973

            # Create calculator
            calc = Calculator(
                get_method("GFN2-xTB"), atomic_numbers, coords_bohr.flatten()
            )

            # Run single point calculation
            res = calc.singlepoint()

            # Get energy in Hartree
            energy_hartree = res.get_energy()

            return {
                "success": True,
                "energy_hartree": energy_hartree,
                "nbasis": "N/A",
                "error_message": "",
            }

        except Exception as e:
            return {"success": False, "error_message": f"{type(e).__name__}: {str(e)}"}
