# Refactoring nat_zacros.py into Separate Repository

## Selected Approach

**Create a new separate repository** under your GitHub account (dauerba/nat_zacros) with the intention to transfer ownership to akandra later for primary maintenance and institutional sharing.

**Directory Structure (After Completion):**
```
c:\Users\a-DJA\GIT\
├── O_Pt111\                          # Original simulation project
│   ├── nat_zacros.py                 # Will be removed after refactoring
│   ├── zacros_functions.py
│   ├── submit_zacros_jobs.py
│   └── ...
│
└── nat_zacros\                       # New separate repository
    ├── nat_zacros\                   # Python package
    │   ├── __init__.py
    │   ├── lattice.py
    │   ├── state.py
    │   ├── trajectory.py
    │   └── rdf.py
    ├── tests\                        # (Future: unit tests)
    ├── examples\                     # (Future: example notebooks)
    ├── setup.py
    ├── README.md
    ├── LICENSE
    └── .gitignore
```

---

## Step-by-Step Refactoring Instructions

### Phase 1: Backup and Setup (5 minutes)

**Step 1: Create backup**
```bash
cd c:\Users\a-DJA\GIT\O_Pt111
cp nat_zacros.py nat_zacros_backup.py
```

**Step 2: Create repository on GitHub**
1. Go to https://github.com/new
2. Repository name: `nat_zacros`
3. Description: "Analysis package for Zacros kinetic Monte Carlo simulations"
4. Public repository (recommended for future sharing)
5. ✅ Add README file
6. ✅ Add .gitignore: Python template
7. ✅ Add license: MIT or Apache-2.0 (recommend MIT for academic tools)
8. Click "Create repository"

**Step 3: Clone locally**
```bash
cd c:\Users\a-DJA\GIT
git clone https://github.com/dauerba/nat_zacros.git
cd nat_zacros
```

**Status after Step 3:** 
- `nat_zacros` directory exists and contains: README.md, .gitignore, LICENSE
- Git repository is initialized and tracking these files

---

### Phase 2: Create Package Structure (45 minutes)

**Step 4: Create package directory**
```bash
# Inside c:\Users\a-DJA\GIT\nat_zacros
mkdir nat_zacros
```

**Step 5: Create `nat_zacros/__init__.py`**

Create file: `c:\Users\a-DJA\GIT\nat_zacros\nat_zacros\__init__.py`

Copy from O_Pt111/nat_zacros.py:
- Lines 1-73: Module docstring (the entire performance guide)

Add after docstring:
```python
"""
[...module docstring from lines 1-73...]
"""

# Import all public classes and functions
from .lattice import lattice
from .state import state
from .trajectory import trajectory
from .rdf import (
    compute_rdf_parallel,
    compute_rdf_parallel_states,
    load_trajectories_parallel
)

# Define what gets imported with "from nat_zacros import *"
__all__ = [
    'lattice',
    'state',
    'trajectory',
    'compute_rdf_parallel',
    'compute_rdf_parallel_states',
    'load_trajectories_parallel'
]

# Package metadata
__version__ = '1.0.0'
__author__ = 'akandra, dauerba'
```

**Step 6: Create `nat_zacros/lattice.py`**

Create file: `c:\Users\a-DJA\GIT\nat_zacros\nat_zacros\lattice.py`

Structure:
```python
"""
Lattice class for FCC(111) surface geometry.
"""

import numpy as np
from pathlib import Path


class lattice:
    """
    [Copy docstring from original if it exists]
    """
    
    # Copy lines 94-383 from O_Pt111/nat_zacros.py
    # This includes:
    # - __init__ method
    # - get_lattice method
    # - __len__ method
    # - apply_pbc method
    # - minimum_image_distance method
    # - pairwise_distances_pbc method
    # - get_nn_distance method
    # - get_cell_area method
    # - __repr__ method
```

**Step 7: Create `nat_zacros/state.py`**

Create file: `c:\Users\a-DJA\GIT\nat_zacros\nat_zacros\state.py`

Structure:
```python
"""
State class for adsorbate configurations on lattice.
"""

import numpy as np
from pathlib import Path


class state:
    """
    [Copy docstring from original if it exists]
    """
    
    # Copy lines 390-519 from O_Pt111/nat_zacros.py
    # This includes:
    # - __init__ method
    # - get_state method
    # - get_coverage method
    # - get_occupied_sites method
    # - get_empty_sites method
    # - get_occupied_coords method
    # - n_ads method
    # - __repr__ method
```

**Step 8: Create `nat_zacros/trajectory.py`**

Create file: `c:\Users\a-DJA\GIT\nat_zacros\nat_zacros\trajectory.py`

Structure:
```python
"""
Trajectory class for sequences of lattice states over time.
"""

import numpy as np
from pathlib import Path
from .state import state  # Note: relative import


class trajectory:
    """
    Container for a sequence of lattice states over time.
    
    [Copy full class docstring from original]
    """
    
    # Copy lines 526-1092 from O_Pt111/nat_zacros.py
    # This includes:
    # - __init__ method
    # - load_trajectory method
    # - add_state method
    # - get_energy_vs_time method
    # - estimate_equilibration method
    # - get_equilibrated_slice method
    # - load_equilibrated_states method
    # - get_g_ref method
    # - get_rdf method
    # - get_cluster_distribution method
    # - get_accessibility_histogram method
    # - get_coverage_vs_time method
    # - __len__ method
    # - __getitem__ method
    # - __repr__ method
```

**Step 9: Create `nat_zacros/rdf.py`**

Create file: `c:\Users\a-DJA\GIT\nat_zacros\nat_zacros\rdf.py`

Structure:
```python
"""
Parallel RDF computation and trajectory loading functions.
"""

import numpy as np
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
from pathlib import Path


# Copy lines 1094-1588 from O_Pt111/nat_zacros.py
# This includes:
# - _compute_single_rdf helper function
# - _compute_state_rdf helper function
# - compute_rdf_parallel function
# - compute_rdf_parallel_states function
# - _load_single_trajectory_equilibrated helper function
# - load_trajectories_parallel function
```

**Step 10: Create `setup.py`**

Create file: `c:\Users\a-DJA\GIT\nat_zacros\setup.py`

```python
from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='nat-zacros',
    version='1.0.0',
    author='akandra, dauerba',
    author_email='your-email@example.com',  # Update with actual email
    description='Analysis package for Zacros kinetic Monte Carlo simulations',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/dauerba/nat_zacros',
    packages=find_packages(),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    python_requires='>=3.8',
    install_requires=[
        'numpy>=1.20',
        'scipy>=1.7',
    ],
    extras_require={
        'parallel': ['tqdm>=4.60'],
        'dev': ['pytest>=6.0', 'pytest-cov>=2.12'],
    },
)
```

**Step 11: Update README.md**

Edit: `c:\Users\a-DJA\GIT\nat_zacros\README.md`

```markdown
# nat_zacros

Python package for analyzing Zacros kinetic Monte Carlo simulations of surface reactions.

## Features

- **Lattice Geometry**: FCC(111) surface lattice with periodic boundary conditions
- **State Management**: Parse and manipulate adsorbate configurations
- **Trajectory Analysis**: Load and analyze KMC trajectories
- **RDF Calculations**: Fast vectorized radial distribution functions with PBC
- **Parallel Processing**: Efficient loading and analysis of multiple trajectories

## Installation

### From GitHub (development)
```bash
pip install git+https://github.com/dauerba/nat_zacros.git
```

### For local development
```bash
git clone https://github.com/dauerba/nat_zacros.git
cd nat_zacros
pip install -e .
```

## Quick Start

```python
from nat_zacros import lattice, trajectory

# Load lattice and trajectory
lat = lattice(dirname='zacros_run_0080')
traj = trajectory(lat, dirname='zacros_run_0080')

# Load trajectory data (equilibrated portion)
traj.load_trajectory(energy_only=True)
traj.load_equilibrated_states(fraction=0.5)

# Compute radial distribution function
r, g = traj.get_rdf(r_max=40.0, dr=0.1)

# Plot
import matplotlib.pyplot as plt
plt.plot(r, g)
plt.xlabel('Distance (Å)')
plt.ylabel('g(r)')
plt.show()
```

## Performance Optimization

The package includes several performance optimizations:

1. **Vectorized distance calculations** (50-100x speedup)
2. **Parallel trajectory loading** (5-10x speedup)
3. **Binary caching with pickle** (100x speedup for repeated analysis)

See module docstring for detailed performance guide.

## Requirements

- Python >= 3.8
- NumPy >= 1.20
- SciPy >= 1.7
- Optional: tqdm (for progress bars in parallel functions)

## Project Structure

```
nat_zacros/
├── nat_zacros/
│   ├── __init__.py       # Package entry point
│   ├── lattice.py        # FCC(111) lattice geometry
│   ├── state.py          # Adsorbate configurations
│   ├── trajectory.py     # State sequences and RDF analysis
│   └── rdf.py            # Parallel computation utilities
├── tests/                # Unit tests (future)
├── examples/             # Example notebooks (future)
├── setup.py              # Package installation
└── README.md
```

## Contributing

This package is part of the O_Pt111 project for studying oxygen adsorption on Pt(111) surfaces.

**Current Maintainers:**
- Primary: akandra (pending repository transfer)
- Developer: dauerba (refactoring and packaging)

## License

MIT License - see LICENSE file for details

## Citation

If you use this package in your research, please cite:
```
[Add citation information when available]
```

## Related Projects

- [Zacros](http://zacros.org/) - Kinetic Monte Carlo software for catalysis
- [O_Pt111](https://github.com/akandra/O_Pt111) - Parent project for Pt(111) simulations
```

---

### Phase 3: Testing (30 minutes)

**Step 12: Commit initial structure**
```bash
cd c:\Users\a-DJA\GIT\nat_zacros
git add .
git commit -m "Initial commit: nat_zacros package structure

Refactored from O_Pt111/nat_zacros.py (1588 lines) into modular package:
- lattice.py: FCC(111) lattice geometry (~300 lines)
- state.py: Adsorbate configurations (~200 lines)
- trajectory.py: State sequences and RDF analysis (~600 lines)
- rdf.py: Parallel computation utilities (~400 lines)

Original development by akandra (O_Pt111 project)
Refactored and packaged by dauerba

Package includes:
- setup.py for pip installation
- README.md with usage examples
- MIT license
- Performance optimization guide in module docstring

Repository created under dauerba account for initial development.
Planned transfer to akandra for primary maintenance and institutional sharing."

git push
```

**Step 13: Test basic import**

Open Python console or new Jupyter notebook:
```python
# Test 1: Basic import
import sys
sys.path.insert(0, r'c:\Users\a-DJA\GIT\nat_zacros')

from nat_zacros import lattice, state, trajectory
print("✓ Import successful")

# Test 2: Create default lattice
lat = lattice()
print(f"✓ Created lattice: {lat}")
print(f"  Sites: {len(lat)}")

# Test 3: Create empty state
st = state(lat)
print(f"✓ Created state: {st}")

# Test 4: Create empty trajectory
traj = trajectory(lat)
print(f"✓ Created trajectory: {traj}")
```

**Step 14: Test with real data**

```python
from pathlib import Path
from nat_zacros import lattice, trajectory

# Use existing simulation data
run_dir = Path(r'c:\Users\a-DJA\GIT\O_Pt111\zacros_run_0080')

# Test loading lattice
lat = lattice(dirname=run_dir)
print(f"✓ Loaded lattice: {len(lat)} sites")
print(f"  Cell area: {lat.get_cell_area():.2f} Å²")

# Test loading trajectory
traj = trajectory(lat, dirname=run_dir)
traj.load_trajectory(energy_only=True)
print(f"✓ Loaded trajectory: {len(traj)} time points")

# Test equilibration
traj.load_equilibrated_states(fraction=0.5)
print(f"✓ Equilibrated states: {len(traj.states)} configurations")

# Test RDF computation
r, g = traj.get_rdf(r_max=40.0, dr=0.1, vectorized=True)
print(f"✓ Computed RDF: {len(r)} distance bins")
print(f"  r range: [{r[0]:.2f}, {r[-1]:.2f}] Å")
```

**Step 15: Test existing notebooks**

1. Open `O_Pt111/test_class.ipynb`
2. Update first cell:
   ```python
   import sys
   sys.path.insert(0, r'c:\Users\a-DJA\GIT\nat_zacros')
   from nat_zacros import lattice, state, trajectory
   ```
3. Run all cells - should work without other changes

**Step 16: Install package locally (optional)**

```bash
cd c:\Users\a-DJA\GIT\nat_zacros
pip install -e .
```

Now you can import from anywhere without sys.path manipulation:
```python
from nat_zacros import lattice, trajectory
```

---

### Phase 4: Update O_Pt111 Repository (15 minutes)

**Step 17: Add nat_zacros as dependency in O_Pt111**

Option A: Link as submodule (keeps repos connected)
```bash
cd c:\Users\a-DJA\GIT\O_Pt111
git submodule add https://github.com/dauerba/nat_zacros.git nat_zacros
git commit -m "Add nat_zacros as submodule"
```

Option B: Use pip install (cleaner, recommended)

Create/update `O_Pt111/requirements.txt`:
```
numpy>=1.20
scipy>=1.7
pandas>=1.3
matplotlib>=3.4
ase>=3.22
networkx>=2.6

# nat_zacros analysis package
nat-zacros @ git+https://github.com/dauerba/nat_zacros.git@main
```

**Step 18: Update copilot-instructions.md**

Edit: `O_Pt111/.github/copilot-instructions.md`

Find the section "Key Components" and update:

```markdown
### Key Components

- **`nat_zacros`** package (separate repository): Core data structures for reading Zacros output
  - Repository: https://github.com/dauerba/nat_zacros
  - `lattice` class: Parses lattice geometry from `lattice_input.dat`/`lattice_output.txt`
  - `state` class: Loads adsorbate configurations from `history_output.txt`
  - `trajectory` class: Manages state sequences and computes RDF with PBC
  - Installation: `pip install git+https://github.com/dauerba/nat_zacros.git`
  
- **`zacros_functions.py`** (1881 lines): Main analysis library with functions for:
  ...
```

Add new section:

```markdown
### Package Migration (2025-12-20)

The `nat_zacros.py` module (1588 lines) was refactored into a separate package:
- **Old location**: `O_Pt111/nat_zacros.py` (removed)
- **New location**: Separate repository at `github.com/dauerba/nat_zacros`
- **Import unchanged**: `from nat_zacros import lattice, state, trajectory`
- **Installation**: Via pip or git submodule (see requirements.txt)

**Rationale for separation:**
- Independent versioning and releases
- Easier sharing with collaborators
- Better maintainability (separate modules vs 1588-line file)
- Planned transfer to akandra for primary maintenance
```

**Step 19: Commit O_Pt111 updates**

```bash
cd c:\Users\a-DJA\GIT\O_Pt111
git add requirements.txt .github/copilot-instructions.md
git commit -m "Update dependencies and docs for nat_zacros package refactoring

- Added nat_zacros package as pip dependency
- Updated copilot-instructions.md with new structure
- Documented package migration rationale

nat_zacros is now maintained as separate repository:
https://github.com/dauerba/nat_zacros"
```

**Step 20: Remove old nat_zacros.py (after confirming everything works)**

```bash
cd c:\Users\a-DJA\GIT\O_Pt111
git rm nat_zacros.py
git commit -m "Remove monolithic nat_zacros.py (replaced by package)

Refactored into separate repository with modular structure.
All functionality preserved with identical API."

git push
```

---

### Phase 5: Future Steps (When akandra is available)

**Step 21: Add akandra as collaborator**
1. Go to https://github.com/dauerba/nat_zacros/settings/access
2. Click "Add people"
3. Search for "akandra"
4. Grant "Write" or "Maintain" access
5. akandra will receive email invitation

**Step 22: Transfer repository ownership (later)**
1. Go to https://github.com/dauerba/nat_zacros/settings
2. Scroll to "Danger Zone"
3. Click "Transfer ownership"
4. Enter new owner: `akandra`
5. akandra must accept transfer

**Step 23: Update URLs after transfer**

In `O_Pt111/requirements.txt`:
```
# Update GitHub URL
nat-zacros @ git+https://github.com/akandra/nat_zacros.git@main
```

---

## Troubleshooting

### Issue: Import fails with ModuleNotFoundError
**Solution:** 
```python
import sys
sys.path.insert(0, r'c:\Users\a-DJA\GIT\nat_zacros')
```
Or install package: `pip install -e c:\Users\a-DJA\GIT\nat_zacros`

### Issue: Relative import error in package files
**Solution:** Use dot notation: `from .state import state`

### Issue: Old nat_zacros.py still loads
**Solution:** Restart Python kernel/console to clear import cache

### Issue: Tests fail in existing notebooks
**Solution:** Update import statement to include package path

---

## Time Estimates

- Phase 1 (Backup & Setup): 5 minutes
- Phase 2 (Create Structure): 45 minutes
- Phase 3 (Testing): 30 minutes
- Phase 4 (Update O_Pt111): 15 minutes
- **Total: ~1.5 hours**

---

## Checklist

- [ ] Step 1: Backup nat_zacros.py
- [ ] Step 2: Create GitHub repository
- [ ] Step 3: Clone locally
- [ ] Step 4: Create package directory
- [ ] Step 5: Create __init__.py
- [ ] Step 6: Create lattice.py
- [ ] Step 7: Create state.py
- [ ] Step 8: Create trajectory.py
- [ ] Step 9: Create rdf.py
- [ ] Step 10: Create setup.py
- [ ] Step 11: Update README.md
- [ ] Step 12: Initial commit
- [ ] Step 13: Test basic import
- [ ] Step 14: Test with real data
- [ ] Step 15: Test existing notebooks
- [ ] Step 16: Install package locally
- [ ] Step 17: Update O_Pt111 dependencies
- [ ] Step 18: Update copilot-instructions.md
- [ ] Step 19: Commit O_Pt111 updates
- [ ] Step 20: Remove old nat_zacros.py
- [ ] Step 21: Add akandra as collaborator
- [ ] Step 22: Transfer repository (future)
- [ ] Step 23: Update URLs after transfer

---

**Last Updated:** 2025-12-20
**Created by:** GitHub Copilot for dauerba
**Purpose:** Guide for refactoring nat_zacros.py into separate package repository