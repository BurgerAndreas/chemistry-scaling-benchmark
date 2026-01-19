# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a computational chemistry benchmarking system that measures time and memory scaling of quantum chemistry methods (UFF, GFN2-xTB, HF, DFT, MP2, CCSD, CCSD(T), FCI) across 85 molecular structures ranging from 10-1000 atoms. The system runs ~2000 calculations with automatic timeout handling, memory tracking, and progressive skipping of expensive calculations.

## Environment Setup

```bash
# Create virtual environment (always use uv)
uv venv

# Install all dependencies
uv pip install pyscf xtb rdkit psutil tqdm numpy typing_extensions
```

## Running Benchmarks

```bash
# Full benchmark (1-2 days runtime, ~2000 calculations)
uv run benchmark_chemistry.py

# Test on single molecule (6 minutes, 32 calculations)
uv run benchmark_chemistry.py --test-mode --molecule xyz_structures/molecule_010atoms.xyz

# Custom timeout (default 600 seconds)
uv run benchmark_chemistry.py --timeout 300

# Custom output file
uv run benchmark_chemistry.py --output my_results.csv
```

## Architecture

### Core Design Pattern: Calculator Abstraction

All quantum chemistry methods inherit from `BaseCalculator` abstract class:

```python
class BaseCalculator(ABC):
    def __init__(self, method_category, method_name, basis_set="N/A")
    @abstractmethod
    def calculate(xyz_file: str) -> Dict[str, Any]  # Must return standardized dict
    def get_method_info() -> Dict[str, str]
```

**Return format contract** (all calculators must follow):
```python
{
    'success': bool,
    'energy_hartree': float,  # Only if success=True
    'nbasis': int or 'N/A',   # Number of basis functions
    'error_message': str      # Only if success=False
}
```

### Three Calculator Implementations

1. **`rdkit_calculator.py`**: UFF force field (basis-free, ~0.01s per molecule)
2. **`xtb_calculator.py`**: GFN2-xTB tight binding (basis-free, ~0.5s per molecule)
3. **`pyscf_calculator.py`**: All basis-set methods (HF, 5 DFT functionals, MP2, CCSD, CCSD(T), FCI)
   - Uses `METHOD_CATEGORIES` dict to map method names to categories
   - Single class handles 10 different quantum chemistry methods via strategy pattern
   - Returns **total energy** (HF + correlation) for MP2/CCSD methods via `e_tot` attribute

### Critical Implementation Details

**PySCF Energy Extraction**:
- MP2/CCSD/CCSD(T): Use `calc.e_tot` (total energy), NOT `calc.kernel()` return value
- `kernel()` returns correlation energy only for MP2/CCSD, which is incorrect for benchmarking
- This is a common pitfall when implementing new wavefunction methods

**RDKit XYZ Parsing**:
- Must call `rdDetermineBonds.DetermineBonds(mol, charge=0)` before `SanitizeMol()`
- XYZ files lack bond information, so connectivity must be determined
- Missing this causes "getNumImplicitHs() called without preceding call to calcImplicitValence()" error

### Timeout + Memory Tracking System

`utils/timeout_runner.py` implements process isolation pattern:

1. Spawns calculation in separate `multiprocessing.Process`
2. Spawns memory monitor process that samples RSS every 100ms via `psutil`
3. Uses `join(timeout=N)` for hard timeout enforcement
4. Returns standardized dict with `success`, `time_seconds`, `peak_memory_mb`, `result`/`error_message`

**Critical**: All calculations run through `run_with_timeout()` wrapper to ensure clean process termination and prevent memory leaks.

### Auto-Skip Logic (BenchmarkRunner)

Two-tier skipping strategy:

1. **Predefined limits** (hard-coded in `should_skip()`):
   - FCI: Skip if natoms > 20 (exponential scaling)
   - CCSD(T): Skip if natoms > 50 (N⁷ scaling)
   - CCSD: Skip if natoms > 100 (N⁶ scaling)

2. **Dynamic timeout tracking** (`timeout_thresholds` dict):
   - Key: `(method_category, method_name, basis_set)` tuple
   - Value: Minimum natoms that caused timeout
   - Updated by `update_timeout_threshold()` after each timeout
   - Prevents redundant expensive calculations

### Output Format

CSV schema (thread-safe writes via `ResultsWriter`):
- `timestamp`, `method_category`, `method_name`, `basis_set`, `molecule_file`, `natoms`
- `success`, `time_seconds`, `peak_memory_mb`, `nbasis`, `energy_hartree`, `error_message`

**Incremental saves**: Each calculation appends to CSV immediately (crash recovery).

## Adding New Methods

To add a new quantum chemistry method:

1. Choose appropriate calculator:
   - Basis-free methods → New calculator class inheriting `BaseCalculator`
   - Basis-set methods → Add to `pyscf_calculator.py` via new `_calculate_X()` method

2. Update `METHOD_CATEGORIES` dict in `pyscf_calculator.py` if using PySCF

3. Add to `run_benchmark()` loop in `benchmark_chemistry.py`:
   ```python
   for xyz_file, natoms in molecules:
       calc = NewCalculator(method_name, basis_set)
       calculations.append((calc, xyz_file, natoms))
   ```

4. **Critical**: Verify energy return value is **total energy**, not correlation energy

5. Add predefined skip limit to `should_skip()` if method scales poorly (e.g., O(N⁶) or worse)

## Expected Scaling Behavior

Use this to validate new methods:

- Force fields: O(N) → linear time, constant memory per atom
- Tight binding: O(N²-N³) → quadratic/cubic time
- HF/DFT: O(N³-N⁴) → cubic/quartic time, memory grows with nbasis²
- MP2: O(N⁵) → quintic scaling
- CCSD: O(N⁶) → timeouts expected at ~100 atoms
- CCSD(T): O(N⁷) → timeouts expected at ~50 atoms
- FCI: O(2^N) → only tractable for ≤20 atoms

## Molecule Dataset

- 85 molecules in `xyz_structures/` directory
- Naming: `molecule_NNNNatoms.xyz` (e.g., `molecule_0100atoms.xyz`)
- Range: 10 atoms (acetone) to 1000 atoms
- Format: Standard XYZ (line 1: natoms, line 2: comment, lines 3+: SYMBOL X Y Z)

## Common Issues

**FCI timeouts on small molecules**: Expected behavior. FCI is exponentially expensive and will timeout even on 10-atom systems with large basis sets.

**MP2/CCSD energies seem wrong**: Check you're using `calc.e_tot` not `calc.kernel()` return value.

**RDKit "Pre-condition Violation"**: Missing `rdDetermineBonds.DetermineBonds()` call before sanitization.

**Memory leaks between calculations**: Ensure all calculations use `run_with_timeout()` wrapper, which spawns isolated processes.
