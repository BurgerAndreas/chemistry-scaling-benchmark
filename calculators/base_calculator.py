"""
Base calculator abstract class.
"""
from abc import ABC, abstractmethod
from typing import Dict, Any
from utils.xyz_parser import parse_xyz


class BaseCalculator(ABC):
    """
    Abstract base class for all calculators.
    """

    def __init__(self, method_category: str, method_name: str, basis_set: str = "N/A"):
        """
        Initialize calculator.

        Args:
            method_category: Category (force_field, tight_binding, HF, DFT, MP2, CCSD, etc.)
            method_name: Specific method name (UFF, GFN2-xTB, RHF, PBE, etc.)
            basis_set: Basis set name (N/A for basis-free methods)
        """
        self.method_category = method_category
        self.method_name = method_name
        self.basis_set = basis_set

    @abstractmethod
    def calculate(self, xyz_file: str) -> Dict[str, Any]:
        """
        Perform calculation on molecule from XYZ file.

        Args:
            xyz_file: Path to XYZ file

        Returns:
            Dictionary with keys:
            - success: bool
            - energy_hartree: float (if success=True)
            - nbasis: int or 'N/A' (number of basis functions)
            - error_message: str (if success=False)
        """
        pass

    def parse_xyz(self, xyz_file: str):
        """
        Parse XYZ file.

        Args:
            xyz_file: Path to XYZ file

        Returns:
            Tuple of (symbols, coordinates, natoms)
        """
        return parse_xyz(xyz_file)

    def get_method_info(self) -> Dict[str, str]:
        """
        Get method information.

        Returns:
            Dictionary with method_category, method_name, basis_set
        """
        return {
            'method_category': self.method_category,
            'method_name': self.method_name,
            'basis_set': self.basis_set
        }
