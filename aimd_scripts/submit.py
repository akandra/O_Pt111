# %%
from pathlib import Path
import sys,os
import numpy as np
import time

p = Path(__file__).resolve().parent

# Get host name
if len(sys.argv)==2:
    host = sys.argv[1]
else:
    print(f"Usage:{sys.argv[0]} [host]")
    exit()
# Check if job script exists
job_script = 'job-'+host+'.sh'
try:
    with open(p / job_script) as f:
        jdata = f.readlines()
except:
    print(f'file {job_script} not found')
    exit()

rng = np.random.default_rng(42)

start_traj = 26
ntrajs = 25
# roll rng
for i in range(0,start_traj-1): rng.random(2)
# loop over trajectories
for i in range(start_traj, start_traj+ntrajs):
    # get random (x,y) coordinates
    x,y = rng.random(2)
    # set initial distance from the surface
    dz = 4 # Angstroem 
    # create a folder for a trajectory
    dir_new = p / str(i)
    try:
        dir_new.mkdir()
    except:
        print(f'directory {dir_new.name} exists')
    else:
        # Copy vasp input files
        for fname in ['POTCAR','KPOINTS','INCAR']:
            try:
                with open(p / fname) as f: data = f.readlines()
                with open(dir_new / fname,'w') as f: f.writelines(data)
            except:
                print(f' {p / fname} not found')
                exit()
        # job script
        try:
            with open(p / job_script) as f: data = f.readlines()
        except:
            print(f' {job_script} not found')
            exit()
        else:
            # change job name
            if host != 'local':
                for il,line in enumerate(data):
                    if line.find('#SBATCH -J')!=-1: data[il] = f'#SBATCH -J {i}unpol\n'
            with open(dir_new / job_script,'w') as f: f.writelines(data)
            
        # Modify POSCAR by changing the projectile's coordinates to (x,y,ztom+dz)
        try:
            with open(p / 'POSCAR.0') as f: data = f.readlines()
        except:
            print(f' POSCAR.0 not found')
            exit()
        else:
            # Get vertical size of the cell
            scale = float(data[1].split()[-1])
            zcell = float(data[4].split()[-1])*scale
            # Get total number of atoms
            nions = sum([int(s) for s in data[6].split()])
            # Get z of the top layer
            ztop = float(data[8+nions-1].split()[2])
            # change the position of the last atom (projectile)
            if data[8].split()[0] == 'Direct': dz = dz/zcell
            data[8+nions] = f'  {x}  {y}  {(ztop + dz):.16f} T T T \n'
            # rewrite POSCAR
            with open(dir_new / 'POSCAR','w') as f: f.writelines(data)

        # Switch to the new folder and run job
        os.chdir(dir_new)
        if host=='local':
            os.system(f'sh {job_script} &')
        else:
            os.system(f'sbatch {job_script}')
        os.system('pwd')
        os.chdir(p)


