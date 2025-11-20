#! /usr/bin/env python-

import sys
import json
import numpy as np
import os
# from itertools import chain
from pathlib import Path
import zacros_cluster_defs as cl

#------------------------------------------------------------------------------
#   Constants class 
#------------------------------------------------------------------------------
class Constants( object ):
	
    def __init__(self):
      self.kb      = 8.6173324E-5     # Boltzmann constant in eV/K
      self.pi      = np.pi

      # Surface specific 
      self.reduced_cell_vectors_pfcc111 = np.array(
                              [ [1.0,          0.0], 
                                [1/2, np.sqrt(3)/2] ] )

      # Cluster catalog data
      self.cluster_dispatch = { '0': cl.cluster_1_site,
                      '1nn': cl.cluster_2_site,
                      '2nn': cl.cluster_3_site_2nn,
                      '3nn': cl.cluster_3_site_3nn,
                      '4nn': cl.cluster_4_site_4nn,
                      '5nn': cl.cluster_4_site_5nn,
                      '6nn': cl.cluster_5_site_6nn,
                      '7nn': cl.cluster_5_site_7nn,
                      '8nn': cl.cluster_5_site_8nn,
                      '9nn': cl.cluster_6_site_9nn,
                    '3body': cl.cluster_3_site_3
                    }

      # 2-body clusters up to 9nn for Pt(111) from Florian's DFT calculations
      self.cl_data_fn_ce = {  '0':  0.0, '1nn':0.123, '2nn':0.029, '3nn':-0.010, '4nn':0.023, 
                           '5nn':0.023, '6nn':0.015, '7nn':0.016, '8nn': 0.030, '9nn':0.014}

      self.cl_data_fn = {  '0':  0.0, '1nn':0.123, '2nn':0.029, '3nn':-0.010, '4nn':0.023, 
                         '5nn':0.023, '6nn':0.015, '7nn':0.016, '8nn': 0.015, '9nn':0.007}

      # 2-body clusters up to 3nn for Pt(111) from Hua's DFT calculations
      self.cl_data_hua = {  '0':  0.0, '1nn':0.191, '2nn':0.034, '3nn':-0.021}

#------------------------------------------------------------------------------
#   Zacros Input Files Class
#------------------------------------------------------------------------------        
class ZacrosInputFiles( object ):
 
    def __init__(self):
        # input file names
        self.simulation_input_file_name       = "simulation_input.dat"
        self.lattice_input_file_name          = "lattice_input.dat"
        self.state_input_file_name            = "state_input.dat"
        self.mechanism_input_file_name        = "mechanism_input.dat"
        self.energetics_input_file_name       = "energetics_input.dat"
        self.path                             = Path("./zacros_jobs")
        self.wdir                             = self.path

        self.run_type                         = ""
        self.n_trajs                          = 1
        self.first_traj                       = 1

        self.header                           = ""
        self.n_ads                            = [0]

        # lattice input parameters
        self.lattice_constant =  1.0   # in Angstroms
        self.reduced_cell_vectors   = np.array([ [1.0, 0.0], [0.0, 1.0] ])
        self.cell_vectors = self.lattice_constant*self.reduced_cell_vectors
        self.repeat_cell = [1,1]   # number of unit cells in directions of lattice vectors
        
        # simulation input parameters
        self.temperature    = 300.0    # in K
        self.snapshots      = 1
        self.max_steps      = 100
        self.max_time       = None
        self.wall_time      = 3600   # in seconds

        # energetics input parameters
        self.cluster_list   = ['0'] # list of cluster types to include
        self.energy_list    = [0.0] # corresponding cluster energies

        self.run_number      = ""       # this is filled in by script
                
    # -------------------------------------------------------------------------
    #  write zacros input files 
    # -------------------------------------------------------------------------
    def write(self, lFile):
        
        # ---------------------------------------------------------------------
        # get the next run number to use in file name
        # ---------------------------------------------------------------------   
        self.run_number = lFile.getNextRunNumber()  
        # leaves the log file open so writing control file and log file are 
        # atomic in a multiprocessor environment. Logfile serves as lock
        
        self.wdir = self.path / f'{self.run_number}'
        try:
          self.wdir.mkdir(parents=True)
        except FileExistsError:
          print('Error: Directory ', self.wdir, ' already exists')
          exit(1)

        # ---------------------------------------------------------------------
        #  write the lattice input file
        # ---------------------------------------------------------------------

        lattice_input_content = [
          f"# {self.header}\n",
          f"lattice periodic_cell\n",
          f"\n",
          f"cell_vectors       # in row format (Angstroms)\n",
          f"\n",
          f"   {self.cell_vectors[0,0]:14.10f}   {self.cell_vectors[0,1]:14.10f} \n",
          f"   {self.cell_vectors[1,0]:14.10f}   {self.cell_vectors[1,1]:14.10f} \n",
          f"\n",
          f"repeat_cell       {self.repeat_cell[0]:.0f} {self.repeat_cell[1]:.0f}\n",
          f"\n",
          f"n_site_types      1\n",
          f"site_type_names   fcc\n",
          f"\n",
          f"n_cell_sites      1\n",
          f"site_types        fcc\n",
          f"\n",
          f"site_coordinates   # fractional coordinates (x,y) in row format\n",
          f"\n",
          f"   0.333333333333333   0.333333333333333\n",
          f"\n",
          f"neighboring_structure\n",
          f"   \n",
          f"   1-1  north\n",
          f"   1-1  east\n",
          f"   1-1  southeast\n",
          f"\n",
          f"end_neighboring_structure\n",
          f"\n",
          f"end_lattice\n",
          ]

        with open(self.wdir / self.lattice_input_file_name, 'w') as f:
            for line in lattice_input_content: f.write(line)
          
        # ---------------------------------------------------------------------
        #  write the simulation input file
        # ---------------------------------------------------------------------

        with open(self.wdir / self.simulation_input_file_name, "w") as f:
          f.write(f"# {self.header}\n")
          f.write(f"\n")
          f.write(f"random_seed               314159265")
          f.write(f"\n")
          f.write(f"temperature               {float(self.temperature)}\n")
          f.write(f"pressure                  1.00\n")
          f.write(f"\n")
          f.write(f"n_gas_species             0\n")
          f.write(f"\n")
          f.write(f"n_surf_species            1\n")
          f.write(f"surf_specs_names          O*\n")
          f.write(f"surf_specs_dent           1\n")
          f.write(f"kmc_propagation_method first_reaction binary_heap\n")
          f.write(f"\n")
          f.write(f"override_array_bounds & & 100 145\n")
          f.write(f"\n")
          f.write(f"snapshots                 on event {self.snapshots}\n")
          f.write(f"process_statistics        off #on realtime 2e-1\n")
          f.write(f"species_numbers           off #on realtime 2e-1\n")
          f.write(f"energetics_lists          off #on realtime 2e-1\n")
          f.write(f"process_lists             off #on realtime 2e-1\n")
          f.write(f"\n")
          f.write(f"event_report              off\n")
          f.write(f"on_sites_seeding_report   off\n")
          f.write(f"\n")
          f.write(f"max_steps                 {self.max_steps}\n")
          if self.max_time is not None:
            f.write(f"max_time                  {self.max_time}\n")
          f.write(f"\n")
          f.write(f"wall_time                 {self.wall_time:.0f} # in seconds\n")
          f.write(f"\n")
          f.write(f"no_restart\n")
          f.write(f"\n")
          f.write(f"# debug_report_processes\n")
          f.write(f"# debug_report_global_energetics\n")
          f.write(f"# debug_check_processes\n")
          f.write(f"# debug_check_lattice\n")
          f.write(f"\n")
          f.write(f"finish\n")

        # ---------------------------------------------------------------------
        #  write the state input file
        # ---------------------------------------------------------------------
        state_input_content = [
            
          f"# {self.header}\n",
          f"\n",
          f"initial_state\n",
          f"\n",
        ]

        line =  [
                  f"seed_multiple O* {self.n_ads[0]}\n",
                  f"  site_types fcc\n",
                  f"end_seed_multiple\n",
                ]
        state_input_content.extend(line)

        state_input_content.append(f"\n")
        state_input_content.append(f"end_initial_state\n")

        with open(self.wdir / self.state_input_file_name, "w") as file:
            for line in state_input_content: file.write(line)

        # ---------------------------------------------------------------------
        #  write mechanism input file
        # ---------------------------------------------------------------------
        mechanism_input_content = [
          f"# {self.header}\n",
          f"\n",
          f"mechanism\n",
          f"\n",
          f"reversible_step O_hopping\n",
          f"  sites 2\n",
          f"  neighboring 1-2\n",
          f"  initial # (entitynumber, species, dentate)\n",
          f"    1 O*    1\n",
          f"	  2 *     1\n",
          f"  final\n",
          f"    2 *     1\n",
          f"	  1 O*    1\n",
          f"  site_types fcc fcc\n",
          f"  pre_expon   5.0e9\n",
          f"  pe_ratio    1.0\n",
          f"  activ_eng   0.43\n",
          f"  prox_factor 0.5\n",
          f"end_reversible_step\n",
          f"\n",
          f"end_mechanism\n",
          f"\n"
        ]

        with open(self.wdir / self.mechanism_input_file_name, "w") as file:
            file.writelines(mechanism_input_content)

        # ---------------------------------------------------------------------
        #  write energetics input file
        # ---------------------------------------------------------------------

        with open(self.wdir / self.energetics_input_file_name, "w") as f:
          f.write('# O at Pt(111)\n')
          f.write('# For structures and values see Dropbox:\n')
          f.write('# "Kinetics of Surface Reactions/zacros/O_Pt111/O_Pt111 structures.pptx"\n')
          f.write('\n')
          f.write('energetics\n')
          f.write('\n')

          for i, s in enumerate(self.cluster_list):
            content = const.cluster_dispatch[s]()
            content.insert(-1, 
              f"  cluster_eng   {self.energy_list[i]:.6f}\n")
            [f.write(line) for line in content]
            f.write('\n')

          f.write('end_energetics\n')

        # ---------------------------------------------------------------------
        #  write the log file and close it
        # ---------------------------------------------------------------------    
        logEntry = [
                    self.run_number, 
                    self.run_type,         
                    self.repeat_cell,
                    self.n_ads,
                    self.temperature,
                    self.cluster_list,
                    self.energy_list
                  ]    
        lFile.writeEntry( logEntry)
          

class logFile ( object):

    #----------------------------------------------------------------------
    # initialize log file: open the file with mode a+
    #    this opens the file if it exists, otherwise creates new file
    #    write a header if opening a new file
    #----------------------------------------------------------------------
    def __init__(self, fname):
      self.fileName = fname
      self.next_run = 1
      header = 'run runType lat_size n_ads temperature clusters energies\n'
      self.log = open(self.fileName, 'a+')
      # if file is at the beginning, it is a new file so add header
      if self.log.tell() == 0:
        self.log.write(header)
      self.log.close()
    
    # -------------------------------------------------------------------------
    #  open log file, get the next run number, leave log file open
    #    the log file is left open so that subesequent operations using it
    #    will be atomic
    # -------------------------------------------------------------------------
    def getNextRunNumber(self):
        # note: this opens the log file and leaves it open
        # file is closed when log file entry is written
        self.log = open(self.fileName, 'a+' )
        self.log.seek(0)
        last_line = self.log.readlines()[-1]

        print()
        print('last line = ', last_line)    
        
        # if last line begins with run it is header --> next run number = 1
        if last_line.split()[0] == 'run':
            return 1
        else:
            return json.loads(last_line)[0] + 1

                              

    # -------------------------------------------------------------------------
    #  write an entry in the log file using json (JavaScript Object Notation)
    #  close the log file
    # -------------------------------------------------------------------------
    def writeEntry(self, entry):
      # Convert any numpy integers to Python integers
      json_entry = []
      for item in entry:
          if isinstance(item, np.integer):
              item = int(item)  # Convert np.int64 to Python int
          elif isinstance(item, list):
              item = [int(x) if isinstance(x, np.integer) else x for x in item]
          json_entry.append(item)

      self.log.write(json.dumps(json_entry) + '\n' )
      self.log.close()
        
        

# -----------------------------------------------------------------------------
#   function to determine if we are running on Linux, OS X or Windows
# -----------------------------------------------------------------------------
def get_platform():
    platforms={
        'linux' : 'Linux',
        'darwin' : 'OS X',
        'win32'  : 'Windows'
        }
    
    if sys.platform not in platforms:
        return sys.platform
    else: 
        return platforms[sys.platform]
          
        
#------------------------------------------------------------------------------
#   Main Program
#------------------------------------------------------------------------------

if "__main__":
    const = Constants()
    print('platform =', get_platform())
    infiles =  ZacrosInputFiles()
    
    # Define simulation conditions
    temperatures = np.arange(100, 201, 100)
    temperature = [70, 200]

    coverages = np.arange(0.1, 0.3, 0.1)
    coverages = [0.02, 0.05, 0.08, 0.1]

    list_of_clusters = [['0'], # ideal lattice gas
                        ['0', '1nn', '2nn'],
                        ['0', '1nn', '2nn', '3nn'],
                        ['0', '1nn', '2nn', '3nn', '4nn', '5nn'],
                        ['0', '1nn', '2nn', '3nn', '4nn', '5nn', '6nn', '7nn', '8nn', '9nn']]
    list_of_energies = [ [ const.cl_data_fn[s] for s in lst ] for lst in list_of_clusters ]

    infiles.path = Path("./zacros_calculations")
    Path(infiles.path).mkdir(parents=True, exist_ok=True)

    infiles.lattice_constant =  2.821135   # in Angstroms
    infiles.reduced_cell_vectors = const.reduced_cell_vectors_pfcc111
    infiles.cell_vectors = infiles.lattice_constant*infiles.reduced_cell_vectors

    infiles.repeat_cell = [40,40]

    infiles.run_type = ""
    infiles.header   = "test"

    # simulation input parameters
    infiles.n_trajs        = 10
    infiles.first_traj     = 1
    infiles.snapshots      = 1000
    infiles.max_steps      = 100000
    infiles.wall_time      = 24*3600  # in seconds
    
    lFile = logFile( infiles.path / "jobs.log" )

    # Looping over conditions and submitting jobs
    for i, cluster_list in enumerate(list_of_clusters):

      infiles.cluster_list = cluster_list
      infiles.energy_list  = list_of_energies[i]

      for T in temperatures:
        for cov in coverages:

          # get temperature and number of adsorbates
          infiles.temperature = T
          infiles.n_ads = [ int(cov*infiles.repeat_cell[0]*infiles.repeat_cell[1]) ]

          # create run directory and write input files
          infiles.write(lFile)

          # write a script to submit jobs as SLURM array
          with open(infiles.wdir / "job.sh", "w") as f:
            f.write("#!/bin/bash\n")
            f.write(f"#SBATCH -J zacros_{infiles.wdir.name}\n")
            f.write(f"#SBATCH --array=1-{infiles.n_trajs}\n")
            f.write(f"#SBATCH --partition=wodtke\n")
            f.write(f"#SBATCH --time=4000:00:00\n")
            f.write(f"#SBATCH --nodes=1\n") 
            f.write(f"#SBATCH --ntasks=1\n")
            f.write(f"# SBATCH --tmp=20G\n")
            f.write(f"# SBATCH --mem=500M\n")
            f.write(f"# SBATCH --mail-user\n")
            f.write(f"# SBATCH --mail-type=END\n")
            f.write(f"# SBATCH -o output\n")
            f.write(f"\n")
            f.write(f"# module purge\n")
            f.write(f"module load shared slurm intel-oneapi compiler-rt ifort tbb mkl mpi\n")
            f.write(f"\n")   
            f.write(f"# export MKL_NUM_THREADS=1                # Use only 1 thread per program\n")
            f.write(f"# export OMP_NUM_THREADS=1                # Use only 1 thread per program\n")

            f.write(f"\n")
            f.write(f"seed=$((314159265 + SLURM_ARRAY_TASK_ID))\n")
            f.write("dir=traj_${SLURM_ARRAY_TASK_ID}\n")
            f.write(f"mkdir -p $dir\n")
            f.write('cp simulation_input.dat $dir/. \n')
            f.write('cp lattice_input.dat $dir/. \n')
            f.write('cp mechanism_input.dat $dir/. \n')
            f.write('cp state_input.dat $dir/. \n')
            f.write('cp energetics_input.dat $dir/. \n')
            f.write('cd $dir \n')
            f.write('sed -i "s/^random_seed.*/random_seed $seed/" "simulation_input.dat"\n')
            f.write('srun /home/akandra/zacros/zacros_4.0/build/zacros.x > "zacros.out"' + '\n')

          # Switch to the working directory
          dir = os.getcwd()
          os.chdir(infiles.wdir)
          #os.system(f'sbatch job.sh')
          os.chdir(dir)

    print('Done submitting jobs')
