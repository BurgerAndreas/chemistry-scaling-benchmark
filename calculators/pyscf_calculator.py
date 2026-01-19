"""
PySCF calculator for HF, DFT, MP2, CCSD, CCSD(T), and FCI methods.
"""

from typing import Dict, Any
from pyscf import gto, scf, dft, mp, cc, fci
from utils.xyz_parser import xyz_to_pyscf_format
from calculators.base_calculator import BaseCalculator


class PySCFCalculator(BaseCalculator):
    """
    Calculator for quantum chemistry methods using PySCF.
    """

    # Map method names to their categories
    METHOD_CATEGORIES = {
        "RHF": "HF",
        "SVWN": "DFT",
        "PBE": "DFT",
        "TPSS": "DFT",
        "B3LYP": "DFT",
        "wB97X": "DFT",
        "MP2": "MP2",
        "CCSD": "CCSD",
        "CCSD(T)": "CCSD(T)",
        "FCI": "FCI",
    }

    def __init__(self, method_name: str, basis_set: str):
        """
        Initialize PySCF calculator.

        Args:
            method_name: Method name (RHF, SVWN, PBE, TPSS, B3LYP, wB97X, MP2, CCSD, CCSD(T), FCI)
            basis_set: Basis set (STO-3G, 6-31G(d), cc-pVDZ)
        """
        if method_name not in self.METHOD_CATEGORIES:
            raise ValueError(f"Unknown method: {method_name}")

        method_category = self.METHOD_CATEGORIES[method_name]

        super().__init__(
            method_category=method_category,
            method_name=method_name,
            basis_set=basis_set,
        )

    def calculate(self, xyz_file: str) -> Dict[str, Any]:
        """
        Perform quantum chemistry calculation.

        Args:
            xyz_file: Path to XYZ file

        Returns:
            Dictionary with success, energy_hartree, nbasis, error_message
        """
        try:
            # Parse XYZ file
            symbols, coords, natoms = self.parse_xyz(xyz_file)

            # Convert to PySCF format
            mol_str = xyz_to_pyscf_format(symbols, coords)

            # Build molecule
            mol = gto.Mole()
            mol.atom = mol_str
            mol.basis = self.basis_set
            mol.unit = "Angstrom"
            mol.spin = 0  # Assume closed-shell
            mol.charge = 0
            mol.verbose = 0  # Suppress output
            mol.build()

            nbasis = mol.nao_nr()

            # Perform calculation based on method
            if self.method_name == "RHF":
                energy = self._calculate_rhf(mol)
            elif self.method_name in ["SVWN", "PBE", "TPSS", "B3LYP", "wB97X"]:
                energy = self._calculate_dft(mol)
            elif self.method_name == "MP2":
                energy = self._calculate_mp2(mol)
            elif self.method_name == "CCSD":
                energy = self._calculate_ccsd(mol)
            elif self.method_name == "CCSD(T)":
                energy = self._calculate_ccsd_t(mol)
            elif self.method_name == "FCI":
                energy = self._calculate_fci(mol)
            else:
                return {
                    "success": False,
                    "error_message": f"Method not implemented: {self.method_name}",
                }

            return {
                "success": True,
                "energy_hartree": energy,
                "nbasis": nbasis,
                "error_message": "",
            }

        except Exception as e:
            return {"success": False, "error_message": f"{type(e).__name__}: {str(e)}"}

    def _calculate_rhf(self, mol) -> float:
        """Calculate RHF energy."""
        mf = scf.RHF(mol)
        mf.conv_tol = 1e-6
        mf.conv_tol_grad = 1e-8
        mf.max_cycle = 100
        energy = mf.kernel()
        if not mf.converged:
            raise RuntimeError("SCF did not converge")
        return energy

    def _calculate_dft(self, mol) -> float:
        """Calculate DFT energy."""
        # Map method names to PySCF functional names
        functional_map = {
            "SVWN": "SVWN",
            "PBE": "PBE",
            "TPSS": "TPSS",
            "B3LYP": "B3LYP",
            "wB97X": "WB97X",
        }

        functional = functional_map[self.method_name]

        mf = dft.RKS(mol)
        mf.xc = functional
        mf.conv_tol = 1e-6
        mf.conv_tol_grad = 1e-8
        mf.max_cycle = 100
        energy = mf.kernel()
        if not mf.converged:
            raise RuntimeError("SCF did not converge")
        return energy

    def _calculate_mp2(self, mol) -> float:
        """Calculate MP2 energy."""
        # First do HF
        mf = scf.RHF(mol)
        mf.conv_tol = 1e-6
        mf.conv_tol_grad = 1e-8
        mf.max_cycle = 100
        mf.kernel()
        if not mf.converged:
            raise RuntimeError("SCF did not converge")

        # Then MP2
        mp2_calc = mp.MP2(mf)
        mp2_calc.kernel()
        # Return total energy (HF + correlation)
        return mp2_calc.e_tot

    def _calculate_ccsd(self, mol) -> float:
        """Calculate CCSD energy."""
        # First do HF
        mf = scf.RHF(mol)
        mf.conv_tol = 1e-6
        mf.conv_tol_grad = 1e-8
        mf.max_cycle = 100
        mf.kernel()
        if not mf.converged:
            raise RuntimeError("SCF did not converge")

        # Then CCSD
        ccsd_calc = cc.CCSD(mf)
        ccsd_calc.conv_tol = 1e-6
        ccsd_calc.max_cycle = 100
        ccsd_calc.kernel()
        if not ccsd_calc.converged:
            raise RuntimeError("CCSD did not converge")
        # Return total energy (HF + correlation)
        return ccsd_calc.e_tot

    def _calculate_ccsd_t(self, mol) -> float:
        """Calculate CCSD(T) energy."""
        # First do HF
        mf = scf.RHF(mol)
        mf.conv_tol = 1e-6
        mf.conv_tol_grad = 1e-8
        mf.max_cycle = 100
        mf.kernel()
        if not mf.converged:
            raise RuntimeError("SCF did not converge")

        # Then CCSD
        ccsd_calc = cc.CCSD(mf)
        ccsd_calc.conv_tol = 1e-6
        ccsd_calc.max_cycle = 100
        ccsd_calc.kernel()
        if not ccsd_calc.converged:
            raise RuntimeError("CCSD did not converge")

        # Calculate (T) correction
        e_t = ccsd_calc.ccsd_t()

        # Total CCSD(T) energy
        energy = ccsd_calc.e_tot + e_t
        return energy

    def _calculate_fci(self, mol) -> float:
        """Calculate FCI energy."""
        # First do HF
        mf = scf.RHF(mol)
        mf.conv_tol = 1e-6
        mf.conv_tol_grad = 1e-8
        mf.max_cycle = 100
        mf.kernel()
        if not mf.converged:
            raise RuntimeError("SCF did not converge")

        # Then FCI
        fci_calc = fci.FCI(mol, mf.mo_coeff)
        energy = fci_calc.kernel()[0]
        return energy
