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


