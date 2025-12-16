# O_Pt111 Project - AI Coding Assistant Guide

## Project Overview
This is a computational materials science project for kinetic Monte Carlo (KMC) simulations of oxygen adsorption on Pt(111) surfaces using Zacros. The codebase analyzes phase diagrams, cluster energetics, and surface reaction kinetics through DFT-parameterized lattice models.

## Core Architecture

### Three Simulation Systems
1. **Zacros KMC** (primary): Kinetic Monte Carlo for surface reactions
   - Input files: `lattice_input.dat`, `mechanism_input.dat`, `energetics_input.dat`, `state_input.dat`
   - Output: `history_output.txt` (trajectory), `lattice_output.txt`, log files
   - Example runs in `zacros_run_0080_upto0nn_0200/`

2. **MMC** (Metropolis Monte Carlo): Grand-canonical ensemble simulations
   - Scripts in `mmc_scripts/`

3. **AIMD** (Ab Initio Molecular Dynamics): DFT-based trajectory analysis
   - Scripts in `aimd_scripts/`

### Key Components
- **`nat_zarcos.py`**: Core data structures (`lattice` and `state` classes) for reading Zacros output
  - `lattice` class: Parses lattice geometry from `lattice_input.dat`/`lattice_output.txt`
  - `state` class: Loads adsorbate configurations from `history_output.txt`
  
- **`zacros_functions.py`** (1881 lines): Main analysis library with functions for:
  - RDF (radial distribution functions)
  - Cluster analysis with PBC (periodic boundary conditions)
  - Visualization using ASE and matplotlib
  - DFT energy extraction and formation energy calculations
  
- **`zacros_cluster_defs.py`**: Generates Zacros cluster definitions for 2-body (up to 9nn) and 3-body interactions

- **`submit_zacros_jobs.py`**: Job submission script with cluster energy data embedded in `Constants` class

## Critical Conventions

### Multi-User Path Management
ALL scripts use a **user-specific path dictionary pattern** for cross-platform/multi-user collaboration:
```python
user_paths = {
    'a-DJA': Path('c:/Users/a-DJA/Dropbox/...'),
    'akandra': Path('/home/akandra/...'),
}
username = os.getenv('USERNAME') if platform.system() == 'Windows' else os.getenv('USER')
userpath = user_paths[username]
```
**Always preserve this pattern when creating new scripts or modifying paths.**

### Surface Science Specifics
- **Lattice**: FCC(111) hexagonal lattice with unit vectors defined in `Constants.reduced_cell_vectors_pfcc111`
- **Site indices**: Zacros uses 1-based indexing; Python conversions use `-1` offset
- **Cluster naming**: `'1nn'`, `'2nn'`, etc. for nearest-neighbor distances; `'1-1-1'`, `'1-2-3'` for 3-body clusters
- **Energy data**: Stored in `Constants` class as dictionaries (`cl_data_fn_ce`, `cl_data_fn_ce_3`) from DFT calculations

### Data Flow Pattern
1. **Setup**: `submit_zacros_jobs.py` generates input files with cluster energetics
2. **Simulation**: External Zacros binary produces `history_output.txt`, log files
3. **Analysis**: 
   - `make_logfile_df.py` → DataFrame from log files
   - `nat_zarcos.py` → Load lattice/state objects
   - `zacros_functions.py` → Compute RDF, cluster distributions, visualizations
4. **Post-processing**: `phase_diagram.py` aggregates results across temperature/chemical potential

## Development Workflows

### Running Simulations
- No automated test suite
- Notebooks (`*.ipynb`) are the primary development/testing environment
- Interactive workflow: create objects → load data → visualize → iterate

### Data Organization
- Calculations stored in user-specific external paths (Dropbox, ownCloud)
- Compressed archives (`.tar.gz`) extracted on-demand in `phase_diagram.py`
- Temporary workspace in `~/kmc_tmp` for analysis

### Jupyter Notebooks
- **`test_class.ipynb`**: Tests `nat_zarcos` classes
- **`phase_diagram.ipynb`**: Phase diagram generation
- **`zacros_analyze.ipynb`**: Trajectory analysis
- **`zacros_cefit.ipynb`**: Cluster expansion fitting
- **`zacros_functions.ipynb`**: Function development/testing

## Dependencies
- **Core**: `numpy`, `pandas`, `matplotlib`, `scipy`
- **Atomistic**: `ase` (Atomic Simulation Environment for visualization)
- **Network**: `networkx` (cluster graph analysis in `clusters_graphs.py`)
- **Optional**: `lmfit` (AIMD fitting), `IPython` (notebook displays)

## Common Pitfalls
- **Path sensitivity**: File existence checks are case-sensitive even on Windows (uses Python, not OS)
- **Index confusion**: Zacros output is 1-indexed; Python arrays are 0-indexed
- **PBC handling**: Unwrap coordinates before computing distances/hulls (see RDF implementation)
- **Large files**: `history_output.txt` can be massive; read selectively using line ranges
- **Data location**: Always check `user_paths` dictionary when scripts fail with path errors

## Key Files for Reference
- [nat_zarcos.py](nat_zarcos.py): Start here for understanding data structures
- [zacros_functions.py](zacros_functions.py): Core analysis functions with docstrings
- [test_class.ipynb](test_class.ipynb): Example usage of `lattice` and `state` classes
- [submit_zacros_jobs.py](submit_zacros_jobs.py): Cluster energy parameters in `Constants` class
