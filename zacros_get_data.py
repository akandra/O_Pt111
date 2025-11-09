#!/usr/bin/env python
#%%
# Imports and initialization of paths and argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import time
import sys

# Script location
p = Path(__file__).absolute().parent

# Check if we are in interactive mode and clear sys.argv[0] if we are to enable use of argparse
if  not ('zacros_get_data.py' in sys.argv[0]):
    sys.argv = ['']

# Argparse setup
parser = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
            description="Scans a directory with Zacros output and extracts necessary data.")

parser.add_argument("source", nargs="?",
                    default="/home/akandra/zacros/kmc_tian/1O",
                    help="absolute path to kMC output")
parser.add_argument("destination", nargs="?",
                    default="kmc_tian/1O", 
                    help="path to write coordinate files relative to script location", )

args = parser.parse_args()

#%%
# Get source and destination folders
in_dir = Path(args.source)
out_dir = p / args.destination

# Get list of trajectories, assuming in_dir refers to zacros
traj_zacros_list = list(in_dir.rglob('history_output.txt'))

if len(traj_zacros_list) > 0: 
    # Zacros
    # Get number of configurations
    nconfs = 0
    with open(traj_zacros_list[0]) as f:
        for line in f:
            if 'configuration' in line:
                nconfs = int(line.split()[1])

    data = np.zeros((nconfs,3))
    start = time.time()
    for itr,tr in enumerate(traj_zacros_list):

        print(f'\rWorking with zacros data...{100*itr/len(traj_zacros_list):.0f}%',end="",flush=True)

        # Get lattice-site coordinates
        lattice = np.loadtxt(tr.parent / 'lattice_output.txt', usecols=(1,2),skiprows=2)
        ncells = len(lattice)

        # Get adsorbate positions
        i = 0
        with open(tr) as f:
            for line in f:
                if 'configuration' in line:
                    data[i,0] = float(line.split()[3])
                    for j in range(ncells):
                        if int(f.readline().split()[2])>0:
                            data[i,1:3] = lattice[j]
                    i += 1

        # Save positions
        output_fname = out_dir / f'x_{tr.parent.name:>06}.dat'
        np.savetxt(output_fname, data, header=f'{'time(s)':10} {'x(Ang)':10} {'y(Ang)':10}')

    end = time.time()
    print(f' done in {end-start:.2f} seconds')

else: 
    # kmc_tian

    # Lattice vectors for fcc(111)
    lvecs = np.array([ [1, 0], [np.cos(np.pi/3), -np.sin(np.pi/3)]])

    traj_tian_list = list(in_dir.glob('*_*.confs'))

    # Get number of configurations
    nconfs = 0
    with open(traj_tian_list[0]) as f:
        for line in f:
            if 'time' in line:
                nconfs += 1

    data = np.zeros((nconfs,3))
    start = time.time()
    for itr,tr in enumerate(traj_tian_list):
        print(f'\rWorking with kmc_tian data...{100*itr/len(traj_tian_list):.0f}%',end="",flush=True)

        # Get adsorbate positions
        i = 0
        with open(tr) as f:
            for line in f:
                if 'time' in line:
                    data[i,0] = float(line.split()[1])
                    f.readline()
                    pos = np.array(f.readline().split()[1:3],dtype=int)
                    data[i,1:3] = np.dot(pos,lvecs)
                    i += 1

        # Save positions
        output_fname = out_dir / f'x_{tr.stem.split('_')[-1]}.dat'
        np.savetxt(output_fname, data, header=f'{'time(s)':10} {'x(Ang)':10} {'y(Ang)':10}')

    end = time.time()
    print(f' done in {end-start:.2f} seconds')


