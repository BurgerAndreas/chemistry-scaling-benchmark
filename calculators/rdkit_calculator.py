"""
RDKit calculator for UFF force field.
"""

from typing import Dict, Any
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds
from calculators.base_calculator import BaseCalculator


class RDKitCalculator(BaseCalculator):
    """
    Calculator for UFF (Universal Force Field) using RDKit.
    """

    def __init__(self):
        super().__init__(
            method_category="force_field", method_name="UFF", basis_set="N/A"
        )

    def calculate(self, xyz_file: str) -> Dict[str, Any]:
        """
        Perform UFF calculation.

        Args:
            xyz_file: Path to XYZ file

        Returns:
            Dictionary with success, energy_hartree, nbasis, error_message
        """
        try:
            # Read XYZ file
            with open(xyz_file, "r") as f:
                xyz_block = f.read()

            # Parse XYZ into RDKit molecule
            mol = Chem.MolFromXYZBlock(xyz_block)

            if mol is None:
                return {"success": False, "error_message": "Failed to parse XYZ file"}

            # Determine connectivity
            rdDetermineBonds.DetermineBonds(mol, charge=0)

            # Sanitize molecule
            Chem.SanitizeMol(mol)

            # Set up UFF force field
            ff = AllChem.UFFGetMoleculeForceField(mol)

            if ff is None:
                return {
                    "success": False,
                    "error_message": "Failed to create UFF force field",
                }

            # Optimize geometry
            converged = ff.Minimize()

            # Get energy in kcal/mol
            energy_kcal = ff.CalcEnergy()

            # Convert to Hartree (1 Hartree = 627.509 kcal/mol)
            energy_hartree = energy_kcal / 627.509

            return {
                "success": True,
                "energy_hartree": energy_hartree,
                "nbasis": "N/A",
                "error_message": "",
            }

        except Exception as e:
            return {"success": False, "error_message": f"{type(e).__name__}: {str(e)}"}
