#! /usr/bin/env python-


"""
Created on Tue Aug 25 13:32:41 2020

@author: dja
"""


#------------------------------------------------------------------------------
#   Constants class 
#------------------------------------------------------------------------------
import sys
import os
from pathlib import Path
import shutil
import subprocess
import numpy as np
import random
import string


class Constants( object ):
	
    def __init__(self):
        self.kb      = 8.6173324E-5     # Boltzmann constant in eV/K
        self.pi      = np.pi

#------------------------------------------------------------------------------
#   Control File Class
#------------------------------------------------------------------------------        
class ControlFile( object ):
 
    def __init__(self):
        # file name components
        self.file_name       = ""
        self.file_name_base  = "test"
        self.run_number      = ""      # this is filled in by script
        self.file_name_type  = "melt" 
        
        # control file parameters
        self.algorithm       = "bkl"			# MC algorithm to use (bkl or mmc)
        self.nlat            = [12,12]  		# size of lattice (n_rows x n_cols)
        self.step_period     = 6			    # step_density = 1/step_period; 
                                                #   0 means no step
        self.adsorbates      = ["O"]            # species on the lattice
        self.coverage        = [0.1]        	# initial partial coverages
        self.temperature     = 1000			    # temperature in K
        self.energy          = "test.energy"    # file name binding and interaction energies
        self.start_conf      = ""               # file name starting configuration
        self.mmc_nsteps      = 10000			# number of MMC steps
        self.mmc_save_period = 10000			# period for output
        self.mmc_hist_period = 0			    # period for histogram  calculation
        self.rdf             = [-1,0.1,100]     # rdf calculation period, bin size and number of bins
        
        self.kmc_ntrajs      = 10
        self.kmc_time        = 0.0001
        self.kmc_nbins       = 1000
        self.kmc_rates       = "test.rates"
        self.kmc_output      = "compressed"
        
        # extra control parameters to control script
        self.temperatures = [600]
        self.n_ads        = [1]
        self.lat_sizes    = [[12,12]]
        self.times        = [ 1e-11,  2e-11,  2e-11,  2e-11,  5e-11,   3e-10   ]
    
                
    # -------------------------------------------------------------------------
    #  write control file 
    # -------------------------------------------------------------------------
    def write(self):
        
        # ---------------------------------------------------------------------
        # construct the file name
        # ---------------------------------------------------------------------   
        self.file_name = self.file_name_base # + '-{:04d}'.format(self.run_number)
        if self.file_name_type != "":
            self.file_name += '-' + self.file_name_type
        
        # ---------------------------------------------------------------------
        #  write the control file
        # ---------------------------------------------------------------------
        with open(self.file_name+'.control', 'w') as f:
            f.write( '! {}\n'.format(self.file_name+'.control') )
            f.write( 'algorithm       {}\n'.format(self.algorithm) )
            f.write( 'nlat            {} {} \n'.format(self.nlat[0], self.nlat[1]) )
            f.write( 'step_period     {}\n'.format(self.step_period) )
            f.write( 'adsorbates      ' + ' '.join(self.adsorbates) + '\n' )
            f.write( 'coverages       ' + ' '.join(['{:.2e}'.format(cov) 
                                                    for cov in self.coverage]) 
                                        +'\n')
            f.write( 'temperature     {}\n'.format(self.temperature))
            f.write( 'energy          {}\n'.format(self.energy))
            f.write( 'rdf             {} {} {} \n'.format(self.rdf[0], self.rdf[1], self.rdf[2]) )
            
            # if self.start_conf is defined and is non null, write start_conf line
            try: 
                if self.start_conf != "":
                    f.write('start_conf      {}\n'.format(self.start_conf))
            except:
                pass   

            if self.algorithm == 'mmc':
                f.write( 'mmc_nsteps      {}\n'.format(self.mmc_nsteps))
                f.write( 'mmc_save_period {}\n'.format(self.mmc_save_period))
                f.write( 'mmc_hist_period {}\n'.format(self.mmc_hist_period))
            
            if self.algorithm == 'bkl':
                f.write( 'kmc_ntrajs      {}\n'.format(self.kmc_ntrajs))
                f.write( 'kmc_time        {}\n'.format(self.kmc_time))
                f.write( 'kmc_nbins       {}\n'.format(self.kmc_nbins))
                f.write( 'kmc_rates       {}\n'.format(self.kmc_rates))
                f.write( 'kmc_output      {}\n'.format(self.kmc_output))


# -----------------------------------------------------------------------------
#   function to generate a random string
# -----------------------------------------------------------------------------
def randStr(chars = string.ascii_uppercase + string.digits, N=10):
	return ''.join(random.choice(chars) for _ in range(N))
        
          
        
#------------------------------------------------------------------------------
#   Main Program
#------------------------------------------------------------------------------

if "__main__":
    const = Constants()
    cFile = ControlFile()
 
    # save script directory so we can copy output files at end of run
    sdir=Path.cwd()
    
    # create a temporary directory
    tdir = Path.home() / "tmp"
    try:
        tdir.mkdir()
    except:
        pass
    
    # create and change directory to a temporary working directory
    wdir = tdir / randStr() 
    try:
        wdir.mkdir()
    except:
        pass   
    os.chdir(wdir)
    
    # copy input files to the working directory
    shutil.copy( sdir / (cFile.file_name_base + ".energy"   ), "." ) 
    shutil.copy( sdir / (cFile.file_name_base + ".rates"), "." )
    
    # loop over lattice sizes (coverages)
    for ilat in range(len(cFile.lat_sizes)):
 
        lsize = cFile.lat_sizes[ilat]
        cFile.nlat = lsize
        print('----------------------------------')
        print('lsize = ', lsize)
        
        for ads in range(len(cFile.n_ads)):
            cFile.coverage[ads] = cFile.n_ads[ads]/(lsize[0]*lsize[1])
        
        # make and change to the execution directory
        exec_dir = Path(str(lsize[0]) +'x'+ str(lsize[1]) + '_' 
                        + str(cFile.temperatures[0]) )
        (wdir / exec_dir).mkdir()
        os.chdir((wdir / exec_dir))
        print('change dir to', (wdir / exec_dir))
        
        #----------------------------------------------------------------------
        #   write and run the melt control file
        #----------------------------------------------------------------------
        cFile.algorithm       = 'mmc'
        cFile.file_name_type  = 'melt'
        cFile.start_conf      = ""
        cFile.temperature     = 1000
        cFile.mmc_nsteps      = 10000
        cFile.mmc_save_period = 10000
        cFile.mmc_hist_period = 0		
        cFile.rdf             = [-1,0.1,100]
        cFile.write()
        print('run melt control file',  cFile.file_name)
#        subprocess.run(['kmc_tian', cFile.file_name])
                
        #----------------------------------------------------------------------
        # loop over temperatures
        #----------------------------------------------------------------------
        for iT in range(len(cFile.temperatures)):
            cFile.temperature = cFile.temperatures[iT]

            if iT == 0:
                # set starting conf filename based on the melt file
                cFile.start_conf = (cFile.file_name +'.confs')
                print('temperature   =', cFile.temperature)
                print('start_conf    =', cFile.start_conf)        
            
            else:
                # set starting conf filename based on the previous run
                cFile.start_conf = ".." / exec_dir / (cFile.file_name +'.confs')
                # shift the execution directory
                exec_dir = wdir / Path(str(lsize[0]) +'x'+ str(lsize[1]) + '_' 
                        + str(cFile.temperature) )
                exec_dir.mkdir()
                os.chdir(exec_dir)
                print('\nchange dir to', exec_dir)
                print('temperature   =', cFile.temperature)
                print('start_conf    =', cFile.start_conf)
        
            #------------------------------------------------------------------
            #   write and run the equilibration control file
            #------------------------------------------------------------------
            cFile.algorithm       = 'mmc'
            cFile.file_name_type  = 'eq'
            cFile.mmc_nsteps      = 1000000
            cFile.mmc_save_period = 10000
            cFile.mmc_hist_period = 1000		
            cFile.rdf             = [1000,0.1,int(cFile.nlat[0]/(2*0.1))]
            cFile.write()
            print('run equilibration control file',  cFile.file_name)
    #        subprocess.run(['kmc_tian', cFile.file_name])

            # -------------------------------------------------------------------------
            # write and run the production control file
            # -------------------------------------------------------------------------
            cFile.algorithm       = 'bkl'
            cFile.start_conf      = ''     # no start configuration for kmc
            cFile.file_name_type  = ''
            cFile.rdf             = [100,0.1,int(cFile.nlat[0]/(2*0.1))]
            cFile.kmc_time        = cFile.times[ilat]
            cFile.write()
            print('run production control file',  cFile.file_name)
#            subprocess.run(['kmc_tian', cFile.file_name ])    
    
    
    
    
    
