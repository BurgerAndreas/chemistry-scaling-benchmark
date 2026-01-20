"""
CSV results writer for benchmark data.
"""

import csv
import os
from datetime import datetime
from typing import Dict, Any
import threading


class ResultsWriter:
    """
    Thread-safe CSV writer for benchmark results.
    """

    def __init__(self, csv_file: str):
        """
        Initialize results writer.

        Args:
            csv_file: Path to CSV output file
        """
        self.csv_file = csv_file
        self.lock = threading.Lock()
        self.fieldnames = [
            "timestamp",
            "method_category",
            "method_name",
            "basis_set",
            "molecule_file",
            "natoms",
            "repetition",
            "success",
            "time_seconds",
            "peak_memory_mb",
            "nbasis",
            "energy_hartree",
            "error_message",
        ]

        # Create file with header if it doesn't exist
        if not os.path.exists(csv_file):
            with open(csv_file, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=self.fieldnames)
                writer.writeheader()

    def write_result(self, result: Dict[str, Any]):
        """
        Write a single result to CSV file.

        Args:
            result: Dictionary with result data
        """
        with self.lock:
            # Ensure all required fields are present
            row = {field: result.get(field, "") for field in self.fieldnames}

            # Add timestamp if not present
            if not row["timestamp"]:
                row["timestamp"] = datetime.now().isoformat()

            with open(self.csv_file, "a", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=self.fieldnames)
                writer.writerow(row)

    def read_results(self) -> list:
        """
        Read all results from CSV file.

        Returns:
            List of result dictionaries
        """
        if not os.path.exists(self.csv_file):
            return []

        with self.lock:
            with open(self.csv_file, "r") as f:
                reader = csv.DictReader(f)
                return list(reader)
