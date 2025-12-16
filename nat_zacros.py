# class definitions for working with Zacros

import numpy as np
from pathlib import Path

#------------------------------------------------------------------------------------------------
# class defining a lattice
#------------------------------------------------------------------------------------------------

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

#------------------------------------------------------------------------------------------------
# class defining a lattice state
#------------------------------------------------------------------------------------------------

class state:
    
    def __init__(self, nsites):
        # default is no species
        self.nsites = nsites
        self.n_gas_species = 0
        self.gas_species_names = []
        self.n_surf_species = 0
        self.surf_species_names = []
        self.surf_species_dent = []

        # Arrays defining the adsorbed species on the lattice
        # with indices corresponding to lattice site indices starting at 1 
        self.ads_ids =    np.zeros(nsites, dtype=int)
        self.occupation = np.zeros(nsites, dtype=int)
        self.dentation =  np.zeros(nsites, dtype=int)


    def get_state(self, dirname, idx=0):

        folder = Path(dirname)

        # Read configuration from history_output.txt file

        try:
            with open(folder / 'history_output.txt', 'r') as f:
                content = f.readlines()    

            for site in range(self.nsites):
                parts = content[7 + idx*(self.nsites+1) + site].split()
                self.ads_ids[site]    = int(parts[1])
                self.occupation[site] = int(parts[2])
                self.dentation[site]  = int(parts[3])
        except:
            print(f'cannot read state_output.txt from {str(folder)}')

    def get_coverage(self):
        """
        Calculate the coverage (fraction of occupied sites).
        
        Returns
        -------
        float
            Fraction of sites that are occupied (0.0 to 1.0)
        """
        return np.count_nonzero(self.occupation) / self.nsites

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

    def get_occupied_coords(self, lattice):
        """
        Get Cartesian coordinates of occupied sites.
        
        Parameters
        ----------
        lattice : lattice object
            Lattice object containing site coordinates
            
        Returns
        -------
        ndarray
            (N, 2) array of coordinates for occupied sites
        """
        mask = self.occupation > 0
        return lattice.coordinates[mask]

    def get_n_adsorbates(self):
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
        n_ads = self.get_n_adsorbates()
        coverage = self.get_coverage()
        return f"state(nsites={self.nsites}, n_adsorbates={n_ads}, coverage={coverage:.3f})"


