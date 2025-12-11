# class definitions for working with Zacros
import numpy as np

# class defining a lattice
class lattice:
    
    def __init__(self):
        # default is an fcc(111) lattice with nearest-neighbor distance of 1.0 Angstrom
        # and a single fcc adsorption site per unit cell
        self.type = "periodic_cell"
        self.r0 = 1.0 # in Angstroms, nearest-neighbor distance
        self.reduced_cell_vectors =  np.array([ [1.0,          0.0], 
                                                [1/2, np.sqrt(3)/2] ])
        self.size = np.array([1, 1])
        self.n_site_types = 1
        self.site_type_names = ['fcc']
        self.n_cell_sites = 1
        self.site_types = ['fcc']
        self.site_coordinates = np.array([[1/3, 1/3]])
        self.neighboring_structure = [ (1,1, 'north'),
                                      (1,1, 'east'),
                                      (1,1, 'southeast') ]
        self.lattice_site_coordinates = np.array([[1/3, 1/3]])
        self.lattice_site_types = np.array([1])
        self.coordination_numbers = np.array([6])
        # list of arrays of nearest neighbors for each lattice site
        self.nn_list = [ np.array([0,0,0,0,0,0]) ] 

# class defining a lattice state
class state:
    
    def __init__(self):
        # default is no gas species and a single adsorbed hydrogen
        self.n_gas_species = 0
        self.gas_species_names = []
        self.n_surf_species = 1
        self.surf_specs_names = ['H*']
        self.surf_specs_dent = [1]
        # Arrays defining the adsorbesd species on the lattice
        # with indices corresponding to lattice site indices starting at 1 )
        self.ads_ids =    np.array([1], dtype=int)
        self.occupation = np.array([1], dtype=int)
        self.dentation =  np.array([1], dtype=int)

