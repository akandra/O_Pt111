#! /usr/bin/env python-


"""
Created on Tue Aug 25 13:32:41 2020

@author: dja
"""


#------------------------------------------------------------------------------
#   Constants class 
#------------------------------------------------------------------------------
import sys
import subprocess
import json
import numpy as np
import time
from itertools import chain
import multiprocessing as mp
from copy import copy, deepcopy

class Constants( object ):
	
    def __init__(self):
        self.kb      = 8.6173324E-5     # Boltzmann constant in eV/K
        self.pi      = np.pi

#------------------------------------------------------------------------------
#   Control File Class
#------------------------------------------------------------------------------        
class mmcControlFile( object ):
 
    def __init__(self):
        # file name components
        self.file_name       = ""
        self.file_name_base  = ""
        self.file_name_type  = "" 
        
        # control file parameters
        self.algorithm       = ""	# MC algorithm to use (bkl or mmc)
        self.nlat            = [0,0]	# size of 2D MC lattice (n_rows x n_cols)
        self.step_period     = 0	# step_density = 1/step_period; 0 means no step
        self.adsorbates      = []	# species on the lattice
        self.coverage        = []	# initial partial coverages per lattice
        self.temperature     = 0	# temperature in K
        self.energy          = ""	# file name with binding and interaction
        self.start_conf      = ""
        self.mmc_nsteps      = 0	# number of MMC steps
        self.mmc_save_period = 0	# period for output
        self.mmc_conf_save   = "false"  # key for saving confs 
        self.mmc_running_avgs_save = "false" # key for saving running averages 
        self.mmc_hist_period = 0	# period for histogram  calculation
        self.rdf             = [-1,0,0] # rdf calculation period, bin size and number of bins
        self.gc_period       = 1        # period for gc-mmc (0 means no gc = canonical mmc)
        self.gc_chempots     = [0]      # chemical potential in eV per species
        self.show_progress   = "false"      # if to show progress
            
        self.run_number      = ""       # this is filled in by script
                
    # -------------------------------------------------------------------------
    #  write control file 
    # -------------------------------------------------------------------------
    def write(self, lFile):
        
        # ---------------------------------------------------------------------
        # get the next run number to use in file name
        # ---------------------------------------------------------------------   
        self.run_number = lFile.getNextRunNumber()  
        # leaves the log file open so writing control file and log file are 
        # atomic in a multiprocessor environment. Logfile serves as lock
        
        # ---------------------------------------------------------------------
        # construct the file name
        # ---------------------------------------------------------------------   
        self.file_name = self.file_name_base + '-{:06d}'.format(self.run_number)
        if self.file_name_type != "":
            self.file_name += '-' + self.file_name_type
        #self.file_name += '.control'
        
        # ---------------------------------------------------------------------
        #  write the control file
        # ---------------------------------------------------------------------
        with open(self.file_name+'.control', 'w') as f:
            f.write( '! {}\n'.format(self.file_name+'.control') )
            f.write( 'algorithm       {}\n'.format(self.algorithm) )
            f.write( 'nlat            {} {} \n'.format(self.nlat[0], self.nlat[1]) )
            f.write( 'step_period     {}\n'.format(self.step_period) )
            f.write( 'adsorbates      ' + ' '.join(self.adsorbates) + '\n' )
            f.write( 'coverages       ' + ' '.join([str(cov) for cov in self.coverage]) +'\n')
            f.write( 'temperature     {}\n'.format(self.temperature))
            f.write( 'energy          {}\n'.format(self.energy))
            
            # if self.start_conf is defined and is non null, write start_conf line
            try: 
                if self.start_conf != "":
                    f.write('start_conf      {}\n'.format(self.start_conf))
            except:
                pass   
            
            f.write( 'mmc_nsteps            {}\n'.format(self.mmc_nsteps))
            f.write( 'mmc_save_period       {}\n'.format(self.mmc_save_period))
            f.write( 'mmc_conf_save         {}\n'.format(self.mmc_conf_save))
            f.write( 'mmc_running_avgs_save {}\n'.format(self.mmc_conf_save))
            f.write( 'mmc_hist_period       {}\n'.format(self.mmc_hist_period))

            f.write( 'rdf             {} {} {} \n'.format(self.rdf[0], self.rdf[1], self.rdf[2]) )
            f.write( 'gc_period       {}\n'.format(self.gc_period))
            f.write( 'gc_chempots     ' + ' '.join([str(cp) for cp in self.gc_chempots]) +'\n')

            f.write( 'show_progress         {}\n'.format(self.show_progress))
            
        # ---------------------------------------------------------------------
        #  write the log file and close it
        # ---------------------------------------------------------------------    
        logEntry = [self.run_number, self.file_name_type, self.nlat,   \
                    self.step_period, self.adsorbates, self.coverage, \
                    self.temperature, self.energy, self.mmc_nsteps,    \
                    self.mmc_save_period, self.mmc_conf_save, \
                    self.mmc_running_avgs_save, \
                    self.mmc_hist_period, self.rdf, \
                    self.gc_period, self.gc_chempots]    
        lFile.writeEntry( logEntry)
          

class logFile ( object):

    #----------------------------------------------------------------------
    # initialize log file: open the file with mode a+
    #    this opens the file if it exists, otherwise creates new file
    #    write a header if opening a new file
    #----------------------------------------------------------------------
    def __init__(self, fnamebase):
        self.fileName = fnamebase + '.log'
        self.next_run = 1
        header = 'run runType nLat nStep ads   cov    Temp     EnergyFile'\
                  + '          Step   Save Confs Ravgs Hist RDF_info gc_period gc_chempots\n'
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
        self.log.write(json.dumps(entry) + '\n' )
        self.log.close()
        
        

# -----------------------------------------------------------------------------
#   function to determine if we are running on Linux, OS X or Windows
# -----------------------------------------------------------------------------
def get_platform():
    platforms={
        'linux1' : 'Linux', 
        'linux2' : 'Linux',
        'darwin' : 'OS X',
        'win32'  : 'Windows'
        }
    
    if sys.platform not in platforms:
        return sys.platform
    else: 
        return platforms[sys.platform]

        
# -----------------------------------------------------------------------------
#   function to do it
# -----------------------------------------------------------------------------
def worker(cFile):
    

    lock.acquire()

    cFile.write(lFile)
    print('running control file',  cFile.file_name)

    lock.release()
    subprocess.run(['./kmc_tian', cFile.file_name])

    return cFile
          
        
#------------------------------------------------------------------------------
#   Main Program
#------------------------------------------------------------------------------

if "__main__":
    const = Constants()
    print('platform =', get_platform())
    cFile_ini =  mmcControlFile()
    lock = mp.Lock()
    
    nprocs = 15
    
    temperatures = range(60, 201, 10)
    coverages = [ [c/1000] for c in range(151, 180, 1) ]

    cFile_ini.nlat            = [36,36]
    cFile_ini.step_period     = 0
    cFile_ini.algorithm       = "mmc"
    cFile_ini.adsorbates      = ["O"]
    cFile_ini.coverage        = [0.01]
    cFile_ini.energy          = "o-pt-fn1.energy"
    cFile_ini.file_name_base  = "fn1"
    cFile_ini.mmc_nsteps      = 1000000
    cFile_ini.mmc_save_period = 1000
    cFile_ini.mmc_hist_period = -1 
    #cFile_ini.mmc_conf_save   = "true"
    #cFile_ini.mmc_running_avgs_save = "true"
    cFile_ini.rdf             = [-1,0.1,50]
    cFile_ini.gc_period       = 0
    #cFile_ini.show_progress   = "true"
    
    lFile = logFile(cFile_ini.file_name_base)
	
    cFiles = []
    
    for iT,T in enumerate(temperatures):
        
        cFiles.append( deepcopy(cFile_ini) )
        cFiles[-1].temperature = T
        cFiles[-1].coverage = coverages[0]

        for cov in coverages[1:]:
            cFiles.append( deepcopy(cFiles[-1]) )
            cFiles[-1].coverage = cov

    mp.Pool(nprocs).map(worker, cFiles)
    
