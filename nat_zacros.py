# class definitions for working with Zacros

import numpy as np
from pathlib import Path

"""
Module: nat_zacros
==================

This module provides classes for working with Zacros simulations:
- `lattice`: FCC(111) surface lattice
- `state`: Adsorbate configuration on the lattice
- `trajectory`: Sequence of states over time

Performance Optimization Guide
-------------------------------
This module includes several performance optimizations for RDF and trajectory analysis.
Understanding when to use each approach is critical for optimal performance.

**Key Optimizations (Recommended for all use cases):**

1. **Vectorized Distance Calculations** (50-100x speedup)
   - Automatically enabled by default in trajectory.get_rdf(vectorized=True)
   - Uses NumPy broadcasting to compute all pairwise distances at once
   - Replaces nested Python loops with compiled NumPy operations
   - No downside, always use this

2. **Binary Caching with pickle** (100x speedup for repeated analysis)
   - Save parsed trajectories to pickle files after first load
   - Subsequent loads read binary instead of parsing text (0.5s vs 60s)
   - Example:
       cache_file = 'trajectories_eq.pkl'
       if not Path(cache_file).exists():
           trajs = [load and parse trajectories]
           with open(cache_file, 'wb') as f:
               pickle.dump(trajs, f)
       else:
           with open(cache_file, 'rb') as f:
               trajs = pickle.load(f)
   - Best for iterative analysis with different parameters

3. **Parallel Loading** (5-10x speedup)
   - Use load_trajectories_parallel() for loading multiple trajectories
   - Each trajectory reads its own file → good I/O parallelism
   - Effective even with optimized parsing

**Advanced Optimizations (Use with caution):**

4. **Parallel RDF Computation** - compute_rdf_parallel(), compute_rdf_parallel_states()
   - WARNING: With vectorization, RDF computation is very fast (~2s for 10 trajectories)
   - Parallelization overhead (process spawn, pickle, IPC) is ~2-4 seconds
   - Only beneficial when computation time >> 20-30 seconds
   - Use cases:
     * compute_rdf_parallel(): >50 trajectories, or very long trajectories
     * compute_rdf_parallel_states(): >100 trajectories with many states each
   - For typical use (10-20 trajectories): sequential computation is FASTER

**Performance Benchmark (typical system: 10 trajectories, ~100 states each, 14 cores):**
- Sequential loading: 64s
- Parallel loading: 6-10s
- RDF computation (sequential, vectorized): 2s
- RDF computation (parallel): 3-5s (slower due to overhead!)

**Recommendation:**
For most users, use vectorized distance calculations and parallel loading.
Only use parallel RDF computation for very large datasets where computation
time significantly exceeds 20-30 seconds.
"""

# --------------------------------------------------------------------------
# ------        Lattice  Class Definition                             ------
# --------------------------------------------------------------------------
#

class lattice:
    
    def __init__(self, dirname=None):

        self.folder = None
        # default is an fcc(111) lattice with nearest-neighbor distance of 1.0 Angstrom
        # and a single fcc adsorption site per unit cell
        self.type = "periodic_cell"
        self.unit_cell_vectors =  np.array([ [1.0,          0.0], 
                                                [1/2, np.sqrt(3)/2] ])
        self.size = np.array([1, 1])
        self.cell_vectors = self.unit_cell_vectors.dot(self.size)
        self.n_site_types = 1
        self.site_type_names = ['fcc']
        self.n_cell_sites = 1
        self.site_types = ['fcc']
        self.fractional_coordinates = np.array([[1/3, 1/3]])
        self.neighboring_structure = [ (1,1, 'north'),
                                      (1,1, 'east'),
                                      (1,1, 'southeast') ]
        self.coordinates = np.array([[1/3, 1/3]])
        self.site_types = np.array([1])
        self.site_coordinations = np.array([6])
        # list of arrays of nearest neighbors for each lattice site
        self.site_nns = [ np.array([0,0,0,0,0,0]) ]

        if dirname is not None:
            self.folder = Path(dirname) 
            self.get_lattice()

    def get_lattice(self):

        if self.folder is None:
            print('nothing to get: lattice folder not defined')
            return

        self.folder = Path(self.folder)
        #
        # Read lattice input file
        #
        try:
            with open(self.folder / 'lattice_input.dat', 'r') as f:
                content = [line for line in f.readlines() if (line.strip() and not line.startswith('#'))]
            content = [line.split('#')[0] for line in content]
            for i,line in enumerate(content):
                ls = line.split()
                if ls[0] == 'lattice':
                    self.type = ls[1]
                if ls[0] == 'cell_vectors':
                    self.unit_cell_vectors = np.array([ [float(x) for x in content[i+1].split()],
                                                        [float(x) for x in content[i+2].split()] ])
                if ls[0] == 'repeat_cell':
                    self.size = np.array([ int(x) for x in ls[1:3] ], dtype=int)
                if ls[0] == 'n_site_types':
                    self.n_site_types = int(ls[1])
                if ls[0] == 'site_type_names':
                    self.site_type_names = ls[1:]
                if ls[0] == 'n_cell_sites':
                    self.n_cell_sites = int(ls[1])
                if ls[0] == 'site_types':
                    self.site_types = ls[1:]
                    # convert to names if given as indices
                    if self.site_types[0].isdigit():
                        self.site_types = [ self.site_type_names[int(x)-1] for x in self.site_types ]
                if ls[0] == 'site_coordinates':
                    self.fractional_coordinates = np.zeros((self.n_cell_sites,2), dtype=float)
                    for j in range(self.n_cell_sites):
                        self.fractional_coordinates[j,:] = np.array([float(x) for x in content[i+1+j].split()[:2]])
                if ls[0] == 'neighboring_structure':
                   self.neighboring_structure = []
                   j = 0
                   while content[i+1+j].split()[0] != 'end_neighboring_structure':
                       parts = content[i+1+j].split()
                       self.neighboring_structure.append( (int(parts[0].split('-')[0]), int(parts[0].split('-')[1]), parts[1]) )
                       j += 1
        except:
            print(f'cannot read lattice_input.dat from {str(self.folder)}')

        #
        # Read lattice output file
        #
        site_coordinates = []
        site_types = []
        site_coordinations = []
        site_nns = []

        try:
            with open(self.folder / 'lattice_output.txt') as f:
                v1 = f.readline().split()[1:3]
                v2 = f.readline().split()[1:3]
                self.cell_vectors = np.array([v1, v2], dtype=float)
                for line in f:
                    ls = line.split()
                    site_coordinates.append(ls[1:3])
                    site_types.append(int(ls[3]))
                    site_coordinations.append(int(ls[4]))
                    # -1 is for python
                    site_nns.append(np.array([ int(ls[5+i])-1 for i in range(int(ls[4]))], dtype=int))

        except:
            print(f'cannot read lattice_output.txt from {str(self.folder)}')

        self.coordinates = np.array(site_coordinates, dtype=float)
        self.site_types = np.array(site_types, dtype=int)
        self.site_coordinations = np.array(site_coordinations, dtype=int)
        self.site_nns = site_nns

    def __len__(self):
        return self.size[0] * self.size[1] * self.n_cell_sites

    def apply_pbc(self, coords):
        """
        Apply periodic boundary conditions to wrap coordinates into primary cell.
        
        Parameters
        ----------
        coords : array_like, shape (N, 2) or (2,)
            Cartesian coordinates to wrap
            
        Returns
        -------
        ndarray
            Wrapped coordinates in primary cell
        """
        coords = np.atleast_2d(coords)
        # Convert to fractional coordinates
        # r_cart = frac @ cell_vectors, so frac = r_cart @ inv(cell_vectors)
        cell_inv = np.linalg.inv(self.cell_vectors)
        frac_coords = coords @ cell_inv
        
        # Wrap to [0, 1)
        frac_coords = frac_coords - np.floor(frac_coords)
        
        # Convert back to Cartesian
        cart_coords = frac_coords @ self.cell_vectors
        
        return cart_coords.squeeze()

    def minimum_image_distance(self, coord1, coord2):
        """
        Calculate minimum image distance between two points with PBC.
        
        Parameters
        ----------
        coord1, coord2 : array_like, shape (2,)
            Cartesian coordinates of two points
            
        Returns
        -------
        float
            Minimum distance respecting periodic boundary conditions
        """
        # Displacement vector
        dr = np.array(coord2) - np.array(coord1)
        
        # Convert to fractional coordinates
        # cell_vectors rows are v1, v2, so: r_cart = frac @ cell_vectors
        # Therefore: frac = r_cart @ inv(cell_vectors)
        cell_inv = np.linalg.inv(self.cell_vectors)
        frac_dr = dr @ cell_inv
        
        # Apply minimum image convention: shift by -1, 0, or +1
        frac_dr = frac_dr - np.rint(frac_dr)
        
        # Convert back to Cartesian
        cart_dr = frac_dr @ self.cell_vectors
        
        return np.linalg.norm(cart_dr)

    def pairwise_distances_pbc(self, coords):
        """
        Calculate all pairwise distances with PBC (vectorized).
        
        Parameters
        ----------
        coords : ndarray, shape (N, 2)
            Cartesian coordinates of N points
            
        Returns
        -------
        distances : ndarray, shape (N, N)
            Matrix of pairwise distances with PBC
            
        Notes
        -----
        This vectorized implementation is much faster than calling 
        minimum_image_distance in nested loops (speedup: ~10-100x for N>100).
        """
        coords = np.asarray(coords)
        n = len(coords)
        
        # Compute all displacement vectors: dr[i,j] = coords[j] - coords[i]
        # Broadcasting: (N, 1, 2) - (1, N, 2) = (N, N, 2)
        dr = coords[np.newaxis, :, :] - coords[:, np.newaxis, :]
        
        # Convert to fractional coordinates
        cell_inv = np.linalg.inv(self.cell_vectors)
        # Reshape for batch matrix multiplication: (N*N, 2) @ (2, 2)
        frac_dr = dr.reshape(-1, 2) @ cell_inv
        
        # Apply minimum image convention
        frac_dr = frac_dr - np.rint(frac_dr)
        
        # Convert back to Cartesian
        cart_dr = frac_dr @ self.cell_vectors
        cart_dr = cart_dr.reshape(n, n, 2)
        
        # Compute norms
        distances = np.linalg.norm(cart_dr, axis=2)
        
        return distances

    def get_nn_distance(self, order=1):
        """
        Get nearest neighbor distance for FCC(111) lattice.
        
        Parameters
        ----------
        order : int
            Neighbor order (1=1nn, 2=2nn, etc.)
            
        Returns
        -------
        float
            Distance to nth nearest neighbor
            
        Notes
        -----
        For FCC(111) with lattice constant a:
        1nn = a, 2nn = sqrt(3)*a, 3nn = 2*a, etc.
        """
        a = np.linalg.norm(self.unit_cell_vectors[0])
        
        # Distance formulas for FCC(111)
        nn_distances = {
            1: a,
            2: np.sqrt(3) * a,
            3: 2 * a,
            4: np.sqrt(7) * a,
            5: 3 * a,
            6: np.sqrt(12) * a,
            7: np.sqrt(13) * a,
            8: 4 * a,
            9: np.sqrt(19) * a,
        }
        
        if order in nn_distances:
            return nn_distances[order]
        else:
            raise ValueError(f"Neighbor order {order} not implemented. Valid orders: 1-9")

    def get_cell_area(self):
        """
        Calculate area of the simulation cell.
        
        Returns
        -------
        float
            Area in square distance units
        """
        # 2D cross product: |v1 × v2| = v1_x * v2_y - v1_y * v2_x
        v1, v2 = self.cell_vectors
        return abs(v1[0] * v2[1] - v1[1] * v2[0])

    def __repr__(self):
        """String representation of lattice"""
        return (f"lattice(type='{self.type}', size={tuple(self.size)}, "
                f"nsites={len(self)}, area={self.get_cell_area():.2f})")





# --------------------------------------------------------------------------
# ------        State Class Definition                                ------
# --------------------------------------------------------------------------

class state:
    

    def __init__(self, lattice, dirname=None):
        self.folder = None

        # Store reference to lattice
        self.lattice = lattice
        
        # default is no species
        self.n_gas_species = 0
        self.gas_species_names = []
        self.n_surf_species = 0
        self.surf_species_names = []
        self.surf_species_dent = []

        # Arrays defining the adsorbed species on the lattice
        # with indices corresponding to lattice site indices starting at 1

        nsites = len(self.lattice) 
        self.ads_ids =    np.zeros(nsites, dtype=int)
        self.occupation = np.zeros(nsites, dtype=int)
        self.dentation =  np.zeros(nsites, dtype=int)

        if dirname is not None:
            self.folder = Path(dirname) 
            self.get_state()


    def get_state(self, idx=0):

        
        # Read configuration from history_output.txt file

        self.folder = Path(self.folder)

        try:
            with open(self.folder / 'history_output.txt', 'r') as f:
                content = f.readlines()    

            nsites = len(self.lattice)
            for site in range(nsites):
                parts = content[7 + idx*(nsites+1) + site].split()
                self.ads_ids[site]    = int(parts[1])
                self.occupation[site] = int(parts[2])
                self.dentation[site]  = int(parts[3])
        except:
            print(f'cannot read history_output.txt from {str(self.folder)}')

    
    def get_coverage(self):
        """
        Calculate the coverage (fraction of occupied sites).
        
        Returns
        -------
        float
            Fraction of sites that are occupied (0.0 to 1.0)
        """
        return np.count_nonzero(self.occupation) / len(self.lattice)

    
    def get_occupied_sites(self):
        """
        Get indices of all occupied sites.
        
        Returns
        -------
        ndarray
            Array of site indices where occupation > 0
        """
        return np.where(self.occupation > 0)[0]

    
    def get_empty_sites(self):
        """
        Get indices of all empty sites.
        
        Returns
        -------
        ndarray
            Array of site indices where occupation == 0
        """
        return np.where(self.occupation == 0)[0]

    
    def get_occupied_coords(self):
        """
        Get Cartesian coordinates of occupied sites.
            
        Returns
        -------
        ndarray
            (N, 2) array of coordinates for occupied sites
        """
        mask = self.occupation > 0
        return self.lattice.coordinates[mask]


    def n_ads(self):
        """
        Get total number of adsorbates on the surface.
        
        Returns
        -------
        int
            Number of occupied sites
        """
        return np.count_nonzero(self.occupation)


    def __repr__(self):
        """String representation of state"""
        coverage = self.get_coverage()
        return f"state(nsites={len(self.lattice)}, n_adsorbates={self.n_ads()}, coverage={coverage:.3f})"


#------------------------------------------------------------------------------------------------
# class defining a trajectory (sequence of states)
#------------------------------------------------------------------------------------------------

class trajectory:
    """
    Container for a sequence of lattice states over time.
    
    Attributes
    ----------
    lattice : lattice object
        The underlying surface lattice
    states : list of state objects
        Sequence of configurations
    times : ndarray
        Time points for each state
    energies : ndarray
        Total energy for each state
    folder : Path
        Directory containing trajectory data
    """
    
# --------------------------------------------------------------------------
# to do: remove lattice from the argument list
# 
    def __init__(self, lattice, dirname=None):
        """
        Initialize trajectory with lattice and optional data folder.
        
        Parameters
        ----------
        lattice : lattice object
            The surface lattice for this trajectory
        dirname : str or Path, optional
            Directory containing history_output.txt
        """
        self.lattice = lattice
        self.states = []
        self.times = []
        self.energies = []
        self.folder = Path(dirname) if dirname else None
        
    def load_trajectory(self, dirname=None, start=0, end=None, step=1, load_energy=True, energy_only=False):
        """
        Load states from history_output.txt.
        
        Parameters
        ----------
        dirname : str or Path, optional
            Override folder location
        start : int, default 0
            First state index to load
        end : int, optional
            Last state index to load (None = all)
        step : int, default 1
            Stride for loading states (applies to configuration indices)
        load_energy : bool, default True
            Whether to extract energy values
        energy_only : bool, default False
            If True, only load time and energy without parsing full state configurations.
            Much faster for energy-only analysis. Recommended to use with step > 1
            for equilibration detection.
        """
        folder = Path(dirname) if dirname else self.folder
        
        if folder is None:
            print('Error: folder not specified')
            return
            
        try:
            if energy_only:
                # Fast path: scan file for configuration headers only
                with open(folder / 'history_output.txt', 'r') as f:
                    idx = 0
                    for line in f:
                        if 'configuration' in line:
                            # Apply start/end/step filters
                            if idx < start:
                                idx += 1
                                continue
                            if end is not None and idx >= end:
                                break
                            if (idx - start) % step != 0:
                                idx += 1
                                continue
                                
                            # Parse time and energy from header
                            parts = line.split()
                            time = float(parts[3])
                            energy = float(parts[5]) if load_energy and len(parts) > 5 else 0.0
                            
                            self.times.append(time)
                            self.energies.append(energy)
                            idx += 1
            else:
                # Full path: load complete state configurations
                with open(folder / 'history_output.txt', 'r') as f:
                    content = f.readlines()
                    
                # Format: header lines, then blocks of (nsites+1) lines per state
                nsites = len(self.lattice)
                n_states = (len(content) - 7) // (nsites + 1)
                
                if end is None:
                    end = n_states
                
                for idx in range(start, min(end, n_states), step):
                    # Parse configuration header line
                    header_line = content[7 + idx * (nsites + 1)]
                    
                    if 'configuration' in header_line:
                        parts = header_line.split()
                        time = float(parts[3])
                        energy = float(parts[5]) if load_energy and len(parts) > 5 else 0.0
                    else:
                        time = idx
                        energy = 0.0
                    
                    # Load full state configuration
                    st = state(self.lattice)
                    st.folder = folder
                    st.get_state(idx=idx)
                    self.states.append(st)
                    self.times.append(time)
                    self.energies.append(energy)
                
        except Exception as e:
            print(f'Error loading trajectory from {str(folder)}: {e}')
            
        self.times = np.array(self.times)
        self.energies = np.array(self.energies)
        
    def add_state(self, state, time=None, energy=None):
        """
        Add a state to the trajectory.
        
        Parameters
        ----------
        state : state object
            Configuration to add
        time : float, optional
            Time point for this state
        energy : float, optional
            Total energy for this state
        """
        self.states.append(state)
        self.times.append(time)
        self.energies.append(energy)
        
    def get_energy_vs_time(self):
        """
        Get energy as a function of time.
        
        Returns
        -------
        times : ndarray
            Time points
        energies : ndarray
            Energy at each time point
        """
        return self.times, self.energies
        
    def estimate_equilibration(self, fraction=0.5, method='fraction'):
        """
        Estimate the index where equilibration begins.
        
        Parameters
        ----------
        fraction : float, default 0.5
            Fraction of trajectory to skip (for method='fraction')
        method : str, default 'fraction'
            Method to use:
            - 'fraction': Skip first fraction of trajectory
            - 'energy_plateau': Detect when energy variance stabilizes (future)
            
        Returns
        -------
        int
            Index where equilibration is considered to start
            
        Notes
        -----
        The default 'fraction' method is simple but effective for most cases.
        Energy typically equilibrates in the first 30-50% of a well-run simulation.
        """
        if method == 'fraction':
            return int(fraction * len(self))
        else:
            raise NotImplementedError(f"Method '{method}' not yet implemented")
            
    def get_equilibrated_slice(self, fraction=0.5, method='fraction'):
        """
        Return a new trajectory containing only equilibrated states.
        
        Parameters
        ----------
        fraction : float, default 0.5
            Fraction of trajectory to skip (passed to estimate_equilibration)
        method : str, default 'fraction'
            Method for equilibration detection
            
        Returns
        -------
        trajectory
            New trajectory object with equilibrated states only
            
        Examples
        --------
        >>> traj = trajectory(lat, dirname)
        >>> traj.load_trajectory()
        >>> traj_eq = traj.get_equilibrated_slice(fraction=0.5)
        >>> r, g = traj_eq.get_rdf()  # RDF only from equilibrated data
        """
        eq_idx = self.estimate_equilibration(fraction=fraction, method=method)
        
        # Create new trajectory with sliced data
        traj_eq = trajectory(self.lattice, self.folder)
        traj_eq.states = self.states[eq_idx:]
        traj_eq.times = self.times[eq_idx:]
        traj_eq.energies = self.energies[eq_idx:]
        
        return traj_eq
        
    def load_equilibrated_states(self, fraction=0.5, method='fraction', dirname=None):
        """
        Reload trajectory with full state data only for equilibrated portion.
        
        This is a two-phase loading strategy:
        1. Uses existing times/energies to determine equilibration point
        2. Reloads only equilibrated configurations with full state data
        
        This avoids loading and parsing non-equilibrated states that would be
        discarded anyway, making RDF/cluster analysis much more efficient.
        
        Parameters
        ----------
        fraction : float, default 0.5
            Fraction of trajectory to skip for equilibration
        method : str, default 'fraction'
            Method for equilibration detection
        dirname : str or Path, optional
            Override folder location
            
        Returns
        -------
        None
            Modifies self.states in place, clearing old states and loading
            only equilibrated configurations.
            
        Examples
        --------
        >>> # Phase 1: Fast energy-only loading
        >>> traj = trajectory(lat, dirname)
        >>> traj.load_trajectory(energy_only=True)
        >>> 
        >>> # Phase 2: Reload equilibrated states for analysis
        >>> traj.load_equilibrated_states(fraction=0.5)
        >>> r, g = traj.get_rdf()  # Now works with full state data
        
        Notes
        -----
        Requires that times/energies are already loaded (from energy_only mode).
        Will clear existing states and reload from file.
        """
        if len(self.times) == 0:
            raise RuntimeError("No trajectory data loaded. Run load_trajectory() first.")
            
        # Determine equilibration index
        eq_idx = self.estimate_equilibration(fraction=fraction, method=method)
        
        # Clear existing states
        self.states = []
        
        # Reload only equilibrated portion with full state data
        folder = Path(dirname) if dirname else self.folder
        
        if folder is None:
            raise RuntimeError('Error: folder not specified')
        
        try:
            # Reload with full state parsing, starting from equilibration point
            # Keep existing times/energies, just populate states
            for idx in range(eq_idx, len(self.times)):
                st = state(self.lattice)
                st.folder = folder
                st.get_state(idx=idx)
                self.states.append(st)
                
        except Exception as e:
            print(f'Error loading equilibrated states from {str(folder)}: {e}')
    
    # ==========================================================================
    # RDF (Radial Distribution Function) Analysis Methods
    # ==========================================================================
    #
    # PERFORMANCE NOTES:
    # ------------------
    # The RDF methods have been heavily optimized through vectorization.
    # Key performance characteristics:
    #
    # 1. get_rdf() with vectorized=True (default):
    #    - Uses lattice.pairwise_distances_pbc() for vectorized distance calculations
    #    - 50-100x faster than nested Python loops
    #    - Typical: 2s for 10 trajectories × 100 states
    #
    # 2. Sequential vs Parallel:
    #    - Sequential loop over trajectories: ~2s for typical datasets
    #    - Parallel RDF functions: ~3-5s (slower due to 2-4s overhead!)
    #    - Recommendation: Use simple sequential loop for RDF computation
    #
    # 3. When parallel RDF helps:
    #    - Only beneficial when computation time >> 30 seconds
    #    - Typically: >50 trajectories or very large systems
    #    - Always benchmark before using parallel functions
    #
    # See module docstring for complete performance analysis.
    # ==========================================================================
    
    def get_g_ref(self, r_max=None, dr=0.1):
        """
        Calculate reference RDF for full lattice (all sites, coverage=1).
        
        This computes the number of neighbors in each distance shell,
        used to normalize the RDF such that g(r)=1 for ideal gas.
        
        Parameters
        ----------
        r_max : float, optional
            Maximum distance for RDF
        dr : float, default 0.1
            Bin width in Angstroms
            
        Returns
        -------
        r_bins : ndarray
            Bin centers
        g_ref : ndarray
            Number of neighbors in each shell (integer counts)
        """
        if r_max is None:
            v1 = self.lattice.cell_vectors[0]
            v2 = self.lattice.cell_vectors[1]
            l1 = np.linalg.norm(v1)
            l2 = np.linalg.norm(v2)
            l3 = np.linalg.norm(v1 + v2)
            r_max = min(l1, l2, l3) / 2.0
        
        # Initialize histogram
        n_bins = int(np.ceil(r_max / dr))
        bin_edges = np.linspace(0.0, r_max, n_bins + 1)
        r_bins = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        counts = np.zeros(n_bins, dtype=int)
        
        # Get all lattice site coordinates
        all_coords = self.lattice.coordinates
        n_sites = len(all_coords)
        
        # Calculate all pairwise distances with PBC
        for i in range(n_sites - 1):
            for j in range(i + 1, n_sites):
                dist = self.lattice.minimum_image_distance(
                    all_coords[i], all_coords[j]
                )
                if 0 < dist <= r_max:
                    bin_idx = int(dist / dr)
                    if bin_idx < n_bins:
                        counts[bin_idx] += 1
        
        # Normalize: 2 * counts / n_sites (factor 2 for unordered pairs)
        g_ref = 2.0 * counts / n_sites
        
        return r_bins, g_ref
        
    def get_rdf(self, r_max=None, dr=0.1, g_ref=None, vectorized=True):
        """
        Calculate radial distribution function averaged over trajectory.
        
        Parameters
        ----------
        r_max : float, optional
            Maximum distance for RDF (default: half of cell diagonal)
        dr : float, default 0.1
            Bin width in Angstroms
        g_ref : ndarray, optional
            Reference RDF for normalization (from full lattice at coverage=1).
            If provided, normalizes by number of neighbors in each shell.
        vectorized : bool, default True
            Use vectorized distance calculations for better performance
            
        Returns
        -------
        r_bins : ndarray
            Bin centers  
        g_r : ndarray
            RDF values normalized such that g(r)=1 for ideal gas
            
        Notes
        -----
        RDF is calculated for occupied sites only and averaged over all states.
        Normalization follows zacros_functions.py: divides counts by g_ref 
        (number of neighbors in each shell) and by coverage.
        """
        if r_max is None:
            # Default: half the minimum cell dimension
            v1 = self.lattice.cell_vectors[0]
            v2 = self.lattice.cell_vectors[1]
            l1 = np.linalg.norm(v1)
            l2 = np.linalg.norm(v2)
            l3 = np.linalg.norm(v1 + v2)
            r_max = min(l1, l2, l3) / 2.0
            
        # Initialize histogram
        n_bins = int(np.ceil(r_max / dr))
        bin_edges = np.linspace(0.0, r_max, n_bins + 1)
        r_bins = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        g_r = np.zeros(n_bins)
        
        # Get average coverage
        avg_coverage = np.mean([s.get_coverage() for s in self.states])
        
        # Accumulate over all states
        for st in self.states:
            occupied_coords = st.get_occupied_coords()
            n_occupied = len(occupied_coords)
            
            if n_occupied < 2:
                continue
            
            counts = np.zeros(n_bins, dtype=int)
            
            if vectorized:
                # Vectorized distance calculation (much faster for large n_occupied)
                distances = self.lattice.pairwise_distances_pbc(occupied_coords)
                # Get upper triangle (no diagonal, no double counting)
                mask = np.triu(np.ones(distances.shape, dtype=bool), k=1)
                valid_dists = distances[mask]
                valid_dists = valid_dists[(valid_dists > 0) & (valid_dists <= r_max)]
                
                # Histogram
                counts, _ = np.histogram(valid_dists, bins=bin_edges)
            else:
                # Original nested loop implementation
                for i in range(n_occupied - 1):
                    for j in range(i + 1, n_occupied):
                        dist = self.lattice.minimum_image_distance(
                            occupied_coords[i], occupied_coords[j]
                        )
                        if 0 < dist <= r_max:
                            bin_idx = int(dist / dr)
                            if bin_idx < n_bins:
                                counts[bin_idx] += 1
            
            # Normalize by g_ref if provided (number of neighbors in each shell)
            if g_ref is not None:
                counts_n = np.zeros_like(counts, dtype=float)
                np.divide(counts, g_ref, out=counts_n, where=g_ref!=0)
                g_r += counts_n / n_occupied / avg_coverage
            else:
                g_r += counts / n_occupied / avg_coverage
        
        # Normalize by number of states
        # Factor of 2 for unordered pairs
        if len(self.states) > 0:
            g_r = 2 * g_r / len(self.states)
                    
        return r_bins, g_r
        
    def get_cluster_distribution(self, nn_cutoff=1):
        """
        Calculate cluster size distribution averaged over trajectory.
        
        Parameters
        ----------
        nn_cutoff : int or float
            Nearest neighbor distance cutoff for clustering.
            If int: nth nearest neighbor distance
            If float: explicit distance in Angstroms
            
        Returns
        -------
        cluster_sizes : ndarray
            Unique cluster sizes
        frequencies : ndarray
            Fraction of time each cluster size appears
            
        Notes
        -----
        Uses connected components algorithm with PBC-aware distances.
        """
        if isinstance(nn_cutoff, int):
            cutoff_dist = self.lattice.get_nn_distance(nn_cutoff) * 1.1  # 10% tolerance
        else:
            cutoff_dist = nn_cutoff
            
        all_clusters = []
        
        for st in self.states:
            occupied_sites = st.get_occupied_sites()
            n_occupied = len(occupied_sites)
            
            if n_occupied == 0:
                continue
                
            # Build adjacency matrix
            occupied_coords = st.get_occupied_coords()
            clusters = []
            visited = np.zeros(n_occupied, dtype=bool)
            
            for i in range(n_occupied):
                if visited[i]:
                    continue
                    
                # Start new cluster with BFS
                cluster = []
                queue = [i]
                visited[i] = True
                
                while queue:
                    current = queue.pop(0)
                    cluster.append(current)
                    
                    # Check neighbors
                    for j in range(n_occupied):
                        if not visited[j]:
                            dist = self.lattice.minimum_image_distance(
                                occupied_coords[current], occupied_coords[j]
                            )
                            if dist < cutoff_dist:
                                visited[j] = True
                                queue.append(j)
                                
                clusters.append(len(cluster))
                
            all_clusters.extend(clusters)
            
        # Calculate distribution
        if len(all_clusters) > 0:
            unique_sizes, counts = np.unique(all_clusters, return_counts=True)
            frequencies = counts / counts.sum()
            return unique_sizes, frequencies
        else:
            return np.array([]), np.array([])
        
    def get_accessibility_histogram(self):
        """
        Calculate histogram of site accessibility (number of vacant nearest neighbors).
        
        Returns
        -------
        accessibility : ndarray
            Number of vacant nearest neighbors (0 to max_coordination)
        frequencies : ndarray
            Fraction of occupied sites with each accessibility
            
        Notes
        -----
        Accessibility measures how many nearest neighbor sites are vacant,
        which affects reactivity and diffusion rates.
        """
        all_accessibility = []
        
        for st in self.states:
            occupied_sites = st.get_occupied_sites()
            
            for site_idx in occupied_sites:
                # Get nearest neighbors for this site
                nn_sites = self.lattice.site_nns[site_idx]
                
                # Count vacant neighbors
                vacant_nn = np.sum(st.occupation[nn_sites] == 0)
                all_accessibility.append(vacant_nn)
                
        # Calculate histogram
        if len(all_accessibility) > 0:
            max_coord = np.max(self.lattice.site_coordinations)
            accessibility = np.arange(max_coord + 1)
            counts = np.zeros(max_coord + 1)
            
            for val in all_accessibility:
                counts[val] += 1
                
            frequencies = counts / counts.sum()
            return accessibility, frequencies
        else:
            return np.array([]), np.array([])
        
    def get_coverage_vs_time(self):
        """
        Get coverage as a function of time.
        
        Returns
        -------
        times : ndarray
            Time points
        coverages : ndarray
            Coverage at each time point
        """
        coverages = np.array([s.get_coverage() for s in self.states])
        return self.times, coverages
        
    def __len__(self):
        """
        Number of configurations in trajectory.
        
        Returns number of states if loaded, otherwise number of time points.
        This allows len() to work correctly for energy_only trajectories.
        """
        if len(self.states) > 0:
            return len(self.states)
        else:
            return len(self.times)
        
    def __getitem__(self, idx):
        """
        Access states by index.
        
        Parameters
        ----------
        idx : int or slice
            Index or slice for states
            
        Returns
        -------
        state or list of states
        
        Examples
        --------
        >>> traj[0]          # First state
        >>> traj[-1]         # Last state  
        >>> traj[10:20:2]    # Every other state from 10 to 20
        """
        return self.states[idx]
        
    def __repr__(self):
        """String representation of trajectory"""
        if len(self) > 0:
            t_range = f"t=[{self.times[0]:.2f}, {self.times[-1]:.2f}]"
        else:
            t_range = "empty"
        return f"trajectory(nstates={len(self)}, {t_range}, lattice={len(self.lattice)} sites)"


# ==============================================================================
# Parallel RDF computation functions
# ==============================================================================
#
# ⚠️  PERFORMANCE WARNING - READ BEFORE USING
# --------------------------------------------
# These parallel functions are rarely beneficial for typical datasets.
# 
# With vectorized distance calculations, RDF computation is very fast (~2s).
# Parallelization overhead (process spawn, pickling, IPC) is ~2-4 seconds.
# 
# RESULT: Sequential computation is FASTER for most use cases.
#
# Use these functions ONLY if:
#   - You have >50 trajectories, OR
#   - Sequential computation takes >30 seconds, AND
#   - You have benchmarked and confirmed speedup
#
# For typical use (10-20 trajectories):
#   Sequential: ~2s
#   Parallel: ~3-5s (SLOWER!)
#
# Recommendation: Use simple sequential loop over trajectories instead.
# ==============================================================================

def _compute_single_rdf(args):
    """
    Helper function for parallel RDF computation (trajectory-level).
    
    Parameters
    ----------
    args : tuple
        (trajectory, r_max, dr, g_ref) tuple for computing RDF
        
    Returns
    -------
    g : ndarray
        RDF values for this trajectory
        
    Notes
    -----
    This is a module-level function to ensure it's pickle-able for multiprocessing.
    """
    traj, r_max, dr, g_ref = args
    r, g = traj.get_rdf(r_max=r_max, dr=dr, g_ref=g_ref)
    return g


def _compute_state_rdf(args):
    """
    Helper function for parallel RDF computation (state-level).
    
    Parameters
    ----------
    args : tuple
        (state, lattice, r_max, dr, g_ref, n_bins) tuple
        
    Returns
    -------
    tuple
        (counts, n_occupied, coverage) for this state
        
    Notes
    -----
    This processes a single state and returns histogram counts for aggregation.
    Used for fine-grained parallelism when you have many states.
    """
    state, lattice, r_max, dr, g_ref, n_bins = args
    
    occupied_coords = state.get_occupied_coords()
    n_occupied = len(occupied_coords)
    coverage = state.get_coverage()
    
    if n_occupied < 2:
        return np.zeros(n_bins), n_occupied, coverage
    
    # Vectorized distance calculation
    distances = lattice.pairwise_distances_pbc(occupied_coords)
    mask = np.triu(np.ones(distances.shape, dtype=bool), k=1)
    valid_dists = distances[mask]
    valid_dists = valid_dists[(valid_dists > 0) & (valid_dists <= r_max)]
    
    # Histogram
    bin_edges = np.linspace(0.0, r_max, n_bins + 1)
    counts, _ = np.histogram(valid_dists, bins=bin_edges)
    
    return counts, n_occupied, coverage


def compute_rdf_parallel(trajectories, r_max=None, dr=0.1, g_ref=None, n_workers=None):
    """
    Compute RDF averaged over multiple trajectories using parallel processing.
    
    ⚠️  WARNING: With vectorized distance calculations, RDF computation is very fast.
    Parallelization overhead (2-4s) often exceeds computation time for typical datasets.
    Only use this function when sequential computation time >> 20-30 seconds.
    
    **When to use:**
    - >50 trajectories
    - Very long trajectories (>1000 states each)
    - Large lattices (>10,000 sites)
    
    **For typical use (10-20 trajectories):** Sequential computation is FASTER.
    Use trajectory.get_rdf() in a simple loop instead.
    
    Parameters
    ----------
    trajectories : list of trajectory objects
        Trajectories to average over. Each trajectory should have states loaded.
    r_max : float, optional
        Maximum distance for RDF calculation (Angstroms). 
        If None, uses minimum cell dimension / 2.
    dr : float, optional
        Bin width for RDF histogram (default: 0.1 Angstrom)
    g_ref : ndarray, optional
        Reference RDF for normalization (from full lattice).
        If None, no normalization by coordination numbers.
    n_workers : int, optional
        Number of parallel workers. If None, uses all available cores.
        
    Returns
    -------
    r : ndarray
        Distance bin centers (Angstroms)
    g_avg : ndarray
        Average RDF over all trajectories
    g_std : ndarray
        Standard deviation of RDF across trajectories
        
    Examples
    --------
    >>> # For typical datasets, use sequential computation (FASTER):
    >>> rdfs = []
    >>> for traj in trajs:
    >>>     r, g = traj.get_rdf(r_max=40.0, dr=0.1, g_ref=g_ref)
    >>>     rdfs.append(g)
    >>> g_avg = np.mean(rdfs, axis=0)
    >>> 
    >>> # Only for very large datasets (>50 trajectories):
    >>> r, g_avg, g_std = compute_rdf_parallel(trajs, r_max=40.0, dr=0.1, 
    ...                                         g_ref=g_ref, n_workers=4)
    
    Notes
    -----
    - Uses ProcessPoolExecutor for true parallel computation
    - Parallelization overhead: ~2-4 seconds (process spawn, pickling, IPC)
    - With vectorization, RDF for 10 trajectories takes ~2s sequentially
    - Therefore parallel version is SLOWER for typical use cases
    - Benchmark your specific case before using this function
    """
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing as mp
    
    # Try to import tqdm for progress bar
    try:
        from tqdm import tqdm
        use_tqdm = True
    except ImportError:
        use_tqdm = False
    
    if len(trajectories) == 0:
        raise ValueError("No trajectories provided")
    
    # Determine number of workers
    if n_workers is None:
        n_workers = mp.cpu_count()
    
    # Prepare arguments for each trajectory
    args_list = [(traj, r_max, dr, g_ref) for traj in trajectories]
    
    # Parallel computation
    print(f"Computing RDF for {len(trajectories)} trajectories using {n_workers} workers...")
    
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        if use_tqdm:
            # With progress bar
            rdfs = list(tqdm(executor.map(_compute_single_rdf, args_list),
                            total=len(trajectories),
                            desc="Computing RDF",
                            unit="traj"))
        else:
            # Without progress bar (fallback)
            rdfs = list(executor.map(_compute_single_rdf, args_list))
    
    # Get distance axis from first trajectory
    r, _ = trajectories[0].get_rdf(r_max=r_max, dr=dr, g_ref=g_ref)
    
    # Compute statistics across trajectories
    rdfs = np.array(rdfs)
    g_avg = np.mean(rdfs, axis=0)
    g_std = np.std(rdfs, axis=0)
    
    print(f"\nSuccessfully computed RDF averaged over {len(trajectories)} trajectories")
    
    return r, g_avg, g_std


def compute_rdf_parallel_states(trajectories, r_max=None, dr=0.1, g_ref=None, n_workers=None):
    """
    Compute RDF with state-level parallelism (better for few trajectories with many states).
    
    ⚠️  WARNING: With vectorized distance calculations, RDF computation is very fast.
    Parallelization overhead (2-4s per trajectory) often exceeds computation time.
    Only use this function when sequential computation time >> 30 seconds.
    
    **When to use:**
    - >100 trajectories with many states each
    - Extremely large systems where per-state computation is slow
    - Already benchmarked and confirmed benefit
    
    **For typical use (10-20 trajectories, ~100 states each):** Sequential is FASTER.
    The overhead from spawning processes for each state dominates.
    
    Parameters
    ----------
    trajectories : list of trajectory objects
        Trajectories to average over. Each trajectory should have states loaded.
    r_max : float, optional
        Maximum distance for RDF calculation (Angstroms). 
    dr : float, optional
        Bin width for RDF histogram (default: 0.1 Angstrom)
    g_ref : ndarray, optional
        Reference RDF for normalization.
    n_workers : int, optional
        Number of parallel workers. If None, uses all available cores.
        
    Returns
    -------
    r : ndarray
        Distance bin centers (Angstroms)
    g_avg : ndarray
        Average RDF over all trajectories
    g_std : ndarray
        Standard deviation of RDF across trajectories
        
    Notes
    -----
    - Parallelizes over individual states instead of trajectories
    - Higher parallelization overhead than compute_rdf_parallel()
    - With vectorization, sequential computation is very fast (~2s for 10x100 states)
    - Parallel overhead exceeds benefit for typical datasets
    - Always benchmark before using in production code
    """
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing as mp
    
    try:
        from tqdm import tqdm
        use_tqdm = True
    except ImportError:
        use_tqdm = False
    
    if len(trajectories) == 0:
        raise ValueError("No trajectories provided")
    
    # Determine r_max if not provided
    if r_max is None:
        traj0 = trajectories[0]
        v1 = traj0.lattice.cell_vectors[0]
        v2 = traj0.lattice.cell_vectors[1]
        l1 = np.linalg.norm(v1)
        l2 = np.linalg.norm(v2)
        l3 = np.linalg.norm(v1 + v2)
        r_max = min(l1, l2, l3) / 2.0
    
    # Initialize histogram
    n_bins = int(np.ceil(r_max / dr))
    bin_edges = np.linspace(0.0, r_max, n_bins + 1)
    r_bins = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    
    # Determine number of workers
    if n_workers is None:
        n_workers = mp.cpu_count()
    
    # Process each trajectory
    rdfs_per_traj = []
    
    for traj_idx, traj in enumerate(trajectories):
        # Prepare arguments for each state in this trajectory
        args_list = [(state, traj.lattice, r_max, dr, g_ref, n_bins) 
                    for state in traj.states]
        
        total_states = len(args_list)
        if total_states == 0:
            continue
        
        print(f"Processing trajectory {traj_idx+1}/{len(trajectories)} ({total_states} states) with {n_workers} workers...")
        
        # Parallel computation across states
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            if use_tqdm:
                results = list(tqdm(executor.map(_compute_state_rdf, args_list, chunksize=max(1, total_states // (n_workers * 4))),
                                total=total_states,
                                desc=f"Traj {traj_idx+1}",
                                unit="state"))
            else:
                results = list(executor.map(_compute_state_rdf, args_list, chunksize=max(1, total_states // (n_workers * 4))))
        
        # Aggregate results for this trajectory
        g_r = np.zeros(n_bins)
        avg_coverage = np.mean([r[2] for r in results])
        
        for counts, n_occupied, coverage in results:
            if n_occupied < 2:
                continue
            
            # Normalize by g_ref if provided
            if g_ref is not None:
                counts_n = np.zeros_like(counts, dtype=float)
                np.divide(counts, g_ref, out=counts_n, where=g_ref!=0)
                g_r += counts_n / n_occupied / avg_coverage
            else:
                g_r += counts / n_occupied / avg_coverage
        
        # Normalize by number of states (factor of 2 for unordered pairs)
        if len(results) > 0:
            g_r = 2 * g_r / len(results)
        
        rdfs_per_traj.append(g_r)
    
    # Compute statistics across trajectories
    rdfs = np.array(rdfs_per_traj)
    g_avg = np.mean(rdfs, axis=0)
    g_std = np.std(rdfs, axis=0)
    
    print(f"\nSuccessfully computed RDF averaged over {len(trajectories)} trajectories")
    print(f"  Total states processed: {sum(len(t.states) for t in trajectories)}")
    
    return r_bins, g_avg, g_std


# ==============================================================================
# Parallel trajectory loading functions
# ==============================================================================

def _load_single_trajectory_equilibrated(args):
    """
    Helper function for parallel trajectory loading.
    
    Parameters
    ----------
    args : tuple
        (lattice, traj_dir, fraction, method) for loading
        
    Returns
    -------
    trajectory
        Trajectory with equilibrated states loaded
        
    Notes
    -----
    Module-level function for multiprocessing pickling.
    """
    lattice, traj_dir, fraction, method = args
    
    # Create trajectory and do two-phase loading
    traj = trajectory(lattice, traj_dir)
    traj.load_trajectory(load_energy=True, energy_only=True)
    traj.load_equilibrated_states(fraction=fraction, method=method)
    return traj
    

def load_trajectories_parallel(lattice, traj_dirs, fraction=0.5, method='fraction', n_workers=None):
    """
    Load multiple trajectories in parallel with equilibrated states.
    
    ✅ RECOMMENDED: This optimization is effective for all use cases.
    I/O operations (file reading/parsing) are slow enough that parallel loading
    provides consistent speedup (5-10x) without excessive overhead.
    
    Unlike parallel RDF computation, loading is I/O-bound and benefits from
    parallelization even for small numbers of trajectories.
    
    Parameters
    ----------
    lattice : lattice object
        Common lattice for all trajectories
    traj_dirs : list of Path
        List of trajectory directory paths
    fraction : float, default 0.5
        Fraction of trajectory to skip for equilibration
    method : str, default 'fraction'
        Equilibration detection method
    n_workers : int, optional
        Number of parallel workers. If None, uses all available cores.
        
    Returns
    -------
    list of trajectory
        Loaded trajectories with equilibrated states
        
    Examples
    --------
    >>> # Sequential loading (slow - 60s for 10 trajectories)
    >>> trajs = []
    >>> for traj_dir in traj_dirs:
    >>>     traj = trajectory(lat, traj_dir)
    >>>     traj.load_trajectory(energy_only=True)
    >>>     traj.load_equilibrated_states(fraction=0.5)
    >>>     trajs.append(traj)
    >>> 
    >>> # Parallel loading (fast - 6-10s for 10 trajectories)
    >>> trajs = load_trajectories_parallel(lat, traj_dirs, fraction=0.5)
    
    Notes
    -----
    - Typical speedup: 5-10x with 10+ cores
    - Effective for all dataset sizes (I/O-bound operations benefit from parallelism)
    - Each trajectory reads its own file → true parallel I/O
    - Overhead is minimal compared to file parsing time
    - Combine with pickle caching for even better repeated-analysis performance
    """
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing as mp
    
    try:
        from tqdm import tqdm
        use_tqdm = True
    except ImportError:
        use_tqdm = False
    
    if len(traj_dirs) == 0:
        return []
    
    # Determine number of workers
    if n_workers is None:
        n_workers = mp.cpu_count()
    
    # Prepare arguments
    args_list = [(lattice, traj_dir, fraction, method) for traj_dir in traj_dirs]
    
    print(f"Loading {len(traj_dirs)} trajectories in parallel using {n_workers} workers...")
    
    # Parallel loading
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        if use_tqdm:
            trajs = list(tqdm(executor.map(_load_single_trajectory_equilibrated, args_list),
                            total=len(traj_dirs),
                            desc="Loading trajectories",
                            unit="traj"))
        else:
            trajs = list(executor.map(_load_single_trajectory_equilibrated, args_list))
    
    print(f"Successfully loaded {len(trajs)} trajectories")
    if len(trajs) > 0:
        if fraction < 1.0:
            print(f"  Example: {len(trajs[0])} states per trajectory (equilibrated)")
        else:
            print(f"  Example: {len(trajs[0])} states per trajectory")
    
    return trajs
