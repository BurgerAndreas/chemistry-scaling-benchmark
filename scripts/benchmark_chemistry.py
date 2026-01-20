#!/usr/bin/env python3
"""
Computational Chemistry Benchmarking System

Benchmarks various computational chemistry methods (force fields, tight binding,
HF, DFT, MP2, CCSD, CCSD(T), FCI) across molecular structures of varying sizes.
"""

import argparse
import os
import glob
from datetime import datetime
from typing import List, Dict, Any
from tqdm import tqdm

from calculators.rdkit_calculator import RDKitCalculator
from calculators.xtb_calculator import XTBCalculator
from calculators.pyscf_calculator import PySCFCalculator
from utils.timeout_runner import run_with_timeout
from utils.results_writer import ResultsWriter
from utils.xyz_parser import parse_xyz


class BenchmarkRunner:
    """
    Main benchmark orchestrator.
    """

    def __init__(self, output_csv: str, timeout_sec: float = 600.0):
        """
        Initialize benchmark runner.

        Args:
            output_csv: Path to output CSV file
            timeout_sec: Timeout in seconds for each calculation
        """
        self.output_csv = output_csv
        self.timeout_sec = timeout_sec
        self.writer = ResultsWriter(output_csv)

        # Track timeout thresholds for skip logic
        # Key: (method_category, method_name, basis_set)
        # Value: max natoms that timed out
        self.timeout_thresholds = {}

    def should_skip(
        self, method_category: str, method_name: str, basis_set: str, natoms: int
    ) -> bool:
        """
        Determine if calculation should be skipped based on timeout history and predefined limits.

        Args:
            method_category: Method category
            method_name: Method name
            basis_set: Basis set
            natoms: Number of atoms

        Returns:
            True if should skip
        """
        # Predefined limits
        predefined_limits = {
            "FCI": 20,
            "CCSD(T)": 50,
            "CCSD": 100,
        }

        # Check predefined limits
        if method_name in predefined_limits:
            if natoms > predefined_limits[method_name]:
                return True

        # Check dynamic timeout threshold
        key = (method_category, method_name, basis_set)
        if key in self.timeout_thresholds:
            threshold = self.timeout_thresholds[key]
            if natoms > threshold:
                return True

        return False

    def update_timeout_threshold(
        self, method_category: str, method_name: str, basis_set: str, natoms: int
    ):
        """
        Update timeout threshold after a timeout occurs.

        Args:
            method_category: Method category
            method_name: Method name
            basis_set: Basis set
            natoms: Number of atoms that timed out
        """
        key = (method_category, method_name, basis_set)
        if key not in self.timeout_thresholds:
            self.timeout_thresholds[key] = natoms
        else:
            # Keep the minimum (most conservative)
            self.timeout_thresholds[key] = min(self.timeout_thresholds[key], natoms)

    def run_calculation(
        self, calculator, xyz_file: str, molecule_file: str, natoms: int, repetition: int
    ) -> Dict[str, Any]:
        """
        Run a single calculation with timeout.

        Args:
            calculator: Calculator instance
            xyz_file: Path to XYZ file
            molecule_file: Molecule filename (for logging)
            natoms: Number of atoms
            repetition: Repetition number

        Returns:
            Result dictionary
        """
        method_info = calculator.get_method_info()

        # Run with timeout
        timeout_result = run_with_timeout(
            calculator.calculate, self.timeout_sec, xyz_file
        )

        # Build result dictionary
        result = {
            "timestamp": datetime.now().isoformat(),
            "method_category": method_info["method_category"],
            "method_name": method_info["method_name"],
            "basis_set": method_info["basis_set"],
            "molecule_file": molecule_file,
            "natoms": natoms,
            "repetition": repetition,
            "success": timeout_result["success"],
            "time_seconds": timeout_result["time_seconds"],
            "peak_memory_mb": timeout_result["peak_memory_mb"],
            "nbasis": "",
            "energy_hartree": "",
            "error_message": "",
        }

        if timeout_result["success"]:
            calc_result = timeout_result["result"]
            result["success"] = calc_result["success"]
            if calc_result["success"]:
                result["energy_hartree"] = calc_result.get("energy_hartree", "")
                result["nbasis"] = calc_result.get("nbasis", "")
            else:
                result["error_message"] = calc_result.get("error_message", "")
        else:
            result["error_message"] = timeout_result.get(
                "error_message", "UNKNOWN_ERROR"
            )

            # Update timeout threshold if timeout occurred
            if result["error_message"] == "TIMEOUT":
                self.update_timeout_threshold(
                    method_info["method_category"],
                    method_info["method_name"],
                    method_info["basis_set"],
                    natoms,
                )

        return result

    def run_benchmark(self, molecule_files: List[str]):
        """
        Run benchmark on all molecule files and methods.

        Args:
            molecule_files: List of XYZ file paths
        """
        # Sort molecules by size
        molecules = []
        for xyz_file in molecule_files:
            try:
                _, _, natoms = parse_xyz(xyz_file)
                molecules.append((xyz_file, natoms))
            except Exception as e:
                print(f"Warning: Could not parse {xyz_file}: {e}")
                continue

        molecules.sort(key=lambda x: x[1])

        # Define all calculations to run
        calculations = []

        # Force field: UFF
        for xyz_file, natoms in molecules:
            calc = RDKitCalculator()
            calculations.append((calc, xyz_file, natoms))

        # Tight binding: GFN2-xTB
        for xyz_file, natoms in molecules:
            calc = XTBCalculator()
            calculations.append((calc, xyz_file, natoms))

        # HF, DFT, MP2, CCSD, CCSD(T), FCI with multiple basis sets
        basis_sets = ["STO-3G", "6-31G(d)", "cc-pVDZ"]
        methods = [
            "RHF",
            "SVWN",
            "PBE",
            "TPSS",
            "B3LYP",
            "wB97X",
            "MP2",
            "CCSD",
            "CCSD(T)",
            "FCI",
        ]

        for method in methods:
            for basis in basis_sets:
                for xyz_file, natoms in molecules:
                    calc = PySCFCalculator(method, basis)
                    calculations.append((calc, xyz_file, natoms))

        # Run all calculations with progress bar
        num_calculations = len(calculations)
        total_runs = num_calculations * 3
        print(f"Running {total_runs} calculations ({num_calculations} unique calculations, 3 repetitions each)...")
        skipped = 0
        failed = 0
        succeeded = 0

        for calc, xyz_file, natoms in tqdm(calculations, desc="Benchmarking"):
            for repetition in range(1, 4):  # Run 3 times
                method_info = calc.get_method_info()
                molecule_file = os.path.basename(xyz_file)

                # Check if should skip
                if self.should_skip(
                    method_info["method_category"],
                    method_info["method_name"],
                    method_info["basis_set"],
                    natoms,
                ):
                    skipped += 1
                    continue

                # Run calculation
                result = self.run_calculation(
                    calc, xyz_file, molecule_file, natoms, repetition
                )

                # Write result
                self.writer.write_result(result)

                if result["success"]:
                    succeeded += 1
                else:
                    failed += 1

        print(f"\nBenchmark complete!")
        print(f"  Succeeded: {succeeded}")
        print(f"  Failed: {failed}")
        print(f"  Skipped: {skipped}")
        print(f"  Total: {total_runs}")
        print(f"\nResults written to: {self.output_csv}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Benchmark computational chemistry methods"
    )
    parser.add_argument(
        "--molecules",
        type=str,
        default="xyz_structures",
        help="Directory containing XYZ molecule files (default: xyz_structures)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output CSV file (default: results/results_{timeout}.csv)",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=600.0,
        help="Timeout in seconds per calculation (default: 600)",
    )
    parser.add_argument(
        "--test-mode", action="store_true", help="Test mode: run on single molecule"
    )
    parser.add_argument(
        "--molecule", type=str, help="Single molecule file for test mode"
    )

    args = parser.parse_args()

    # Set default output file if not specified
    if args.output is None:
        os.makedirs("results", exist_ok=True)
        args.output = f"results/results_{int(args.timeout)}.csv"

    # Find molecule files
    if args.test_mode and args.molecule:
        molecule_files = [args.molecule]
    else:
        molecule_pattern = os.path.join(args.molecules, "*.xyz")
        molecule_files = glob.glob(molecule_pattern)

        if not molecule_files:
            print(f"Error: No XYZ files found in {args.molecules}")
            return 1

    print(f"Found {len(molecule_files)} molecule files")

    # Run benchmark
    runner = BenchmarkRunner(args.output, args.timeout)
    runner.run_benchmark(molecule_files)

    return 0


if __name__ == "__main__":
    exit(main())
