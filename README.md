# Computational Chemistry Benchmarking System

A comprehensive benchmarking suite for measuring time and memory scaling of computational chemistry methods across molecular structures of varying sizes (10-1000 atoms).

## Features

- **Multiple Method Categories**: Force fields, tight binding, Hartree-Fock, DFT, MP2, CCSD, CCSD(T), and FCI
- **Automatic Timeout Handling**: 10-minute timeout per calculation with clean process termination
- **Memory Tracking**: Real-time peak memory monitoring using psutil
- **Progressive Skipping**: Automatically skips expensive calculations on larger molecules based on timeout history
- **Incremental Results**: Saves results after each calculation to enable crash recovery
- **Comprehensive Coverage**: ~2000 total calculations across 85 molecular structures

## Implemented Methods

### 1. Force Field
- **UFF** (Universal Force Field) - RDKit implementation

### 2. Tight Binding
- **GFN2-xTB** - Semi-empirical tight binding method

### 3. Hartree-Fock
- **RHF** (Restricted Hartree-Fock)
- Basis sets: STO-3G, 6-31G(d), cc-pVDZ

### 4. DFT (Density Functional Theory)
All with 3 basis sets (STO-3G, 6-31G(d), cc-pVDZ):
- **SVWN** - LDA functional
- **PBE** - GGA functional
- **TPSS** - Meta-GGA functional
- **B3LYP** - Hybrid functional
- **ωB97X** - Range-separated functional

### 5. Wavefunction Methods
All with 3 basis sets (STO-3G, 6-31G(d), cc-pVDZ):
- **MP2** - Møller-Plesset 2nd order perturbation theory
- **CCSD** - Coupled Cluster Singles Doubles
- **CCSD(T)** - CCSD with perturbative triples
- **FCI** - Full Configuration Interaction

## Installation

```bash
# Create virtual environment
uv venv

# Install dependencies
uv pip install pyscf xtb rdkit psutil tqdm numpy typing_extensions
```

## Usage

### Full Benchmark Run
```bash
uv run scripts/benchmark_chemistry.py
```

### Custom Options
```bash
# Specify molecule directory
uv run scripts/benchmark_chemistry.py --molecules xyz_structures

# Custom output file
uv run scripts/benchmark_chemistry.py --output my_results.csv

# Adjust timeout (in seconds)
uv run scripts/benchmark_chemistry.py --timeout 600

# Test mode on single molecule
uv run scripts/benchmark_chemistry.py --test-mode --molecule xyz_structures/molecule_010atoms.xyz
```

## Output Format

Results are saved to CSV with the following columns:

- `timestamp`: ISO-8601 timestamp
- `method_category`: force_field | tight_binding | HF | DFT | MP2 | CCSD | CCSD(T) | FCI
- `method_name`: Specific method (UFF, GFN2-xTB, RHF, PBE, etc.)
- `basis_set`: N/A | STO-3G | 6-31G(d) | cc-pVDZ
- `molecule_file`: XYZ filename
- `natoms`: Number of atoms
- `success`: True | False
- `time_seconds`: Wall-clock time (600.0 for timeout)
- `peak_memory_mb`: Peak RSS memory in MB
- `nbasis`: Number of basis functions
- `energy_hartree`: Total energy in Hartree
- `error_message`: Error description if failed

## Architecture

```
/ssd/Code/scaling/
├── scripts/
│   ├── benchmark_chemistry.py      # Main orchestrator
│   └── figure_tradeoff.py          # Generate tradeoff plots
├── calculators/
│   ├── __init__.py
│   ├── base_calculator.py          # Abstract base class
│   ├── rdkit_calculator.py         # UFF force field
│   ├── xtb_calculator.py           # GFN2-xTB
│   └── pyscf_calculator.py         # HF, DFT, MP2, CCSD, CCSD(T), FCI
├── utils/
│   ├── __init__.py
│   ├── timeout_runner.py           # Timeout + memory tracking
│   ├── xyz_parser.py               # Parse XYZ files
│   └── results_writer.py           # Thread-safe CSV writing
├── datagen/
│   ├── fetch_molecules.py          # Fetch molecules from databases
│   └── generate_small_molecules.py # Generate small test molecules
├── xyz_structures/                 # Molecular structure files (.xyz)
└── results/
    ├── results_{timeout}.csv       # Output data
    └── *.png                       # Generated plots
```

## Performance Characteristics

Based on test run with 10-atom molecule:

| Method | Basis Set | Time (s) | Memory (MB) | Status |
|--------|-----------|----------|-------------|--------|
| UFF | N/A | 0.01 | 80 | ✓ |
| GFN2-xTB | N/A | 0.5 | 101 | ✓ |
| RHF | STO-3G | 1.1 | 104 | ✓ |
| RHF | 6-31G(d) | 2.1 | 136 | ✓ |
| RHF | cc-pVDZ | 2.7 | 162 | ✓ |
| B3LYP | STO-3G | 5.1 | 202 | ✓ |
| B3LYP | 6-31G(d) | 7.5 | 346 | ✓ |
| B3LYP | cc-pVDZ | 8.0 | 433 | ✓ |
| MP2 | STO-3G | 1.5 | 118 | ✓ |
| MP2 | cc-pVDZ | 3.1 | 179 | ✓ |
| CCSD | STO-3G | 4.4 | 114 | ✓ |
| CCSD(T) | cc-pVDZ | 13.0 | 554 | ✓ |
| FCI | STO-3G | 60.0 | 12442 | Timeout |

## Auto-Skip Logic

The system implements intelligent skipping:

- **Predefined limits**:
  - FCI: Only molecules ≤20 atoms
  - CCSD(T): Only molecules ≤50 atoms
  - CCSD: Only molecules ≤100 atoms

- **Dynamic timeout tracking**: If a method times out on N atoms, automatically skip all molecules >N atoms

## Expected Scaling

- **Force fields (UFF)**: O(N) - linear scaling
- **Tight binding (GFN2-xTB)**: O(N²-N³) - quadratic to cubic
- **HF/DFT**: O(N³-N⁴) - cubic to quartic
- **MP2**: O(N⁵) - quintic scaling
- **CCSD**: O(N⁶) - hexic scaling
- **CCSD(T)**: O(N⁷) - septic scaling
- **FCI**: O(2^N) - exponential (intractable for >~20 atoms)

## Test Results

Test run on 10-atom molecule (acetone):
- **Total calculations**: 32
- **Succeeded**: 29 (91%)
- **Failed**: 3 (9% - all FCI timeouts)
- **Runtime**: ~6 minutes

All methods except FCI work correctly. FCI timeout is expected behavior.

## Full Benchmark Estimate

For 85 molecules:
- **Estimated calculations**: ~2000 (with auto-skip)
- **Estimated time**: 1-2 days
- **Output size**: ~2000 rows in CSV

## Notes

- Uses spin-restricted calculations (RHF, RKS) assuming closed-shell molecules
- SCF convergence: energy 1e-6, density 1e-8
- Memory tracking samples every 100ms
- Results append to CSV after each calculation for crash recovery
- All energies reported in Hartree atomic units
