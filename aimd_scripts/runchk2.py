#!/usr/bin/env python3
# consider to make trajectory list dynamic using df
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

# Get list of trajectories to oversee
#traj_list = sorted([ d for d in p.iterdir() if d.is_dir()])
traj_list = []
for d in p.iterdir():
    if d.is_dir():
        try:
            int(d.name)
            traj_list.append(d)
        except: pass
# Max number of ionic steps
nsteps_max = 400
# Max time for a run on a cluster
max_time = 22*3600 # seconds
# how often to check the trajectory status
chk_time = 10*60 # seconds
# Maximum allowd existence of STOPCAR
max_stopcar_life = 3*3600 # seconds
# width of the region below z_stop (in Angstrom) 
# where traj is checked for scattering
z_buffer = 1.0
z_stop = 4
# Unit conversion factor: amu*(A/fs)^2 to eV
feV = 103.6382 
# state log file name
runchk_fname = p / "runchk.info"

# Initialize important arrays
start_time   = np.zeros(len(traj_list))
current_time = np.zeros(len(traj_list))
restart_counter = np.zeros(len(traj_list), dtype=int)
z = np.zeros(len(traj_list))
status = [ 'not_started' for _ in range(len(traj_list)) ]
done = np.full(len(traj_list), False)
nsteps = np.zeros(len(traj_list), dtype=int)

# Initial loop to get trajectories status
for i, wdir in enumerate(traj_list):

    # check if OURCAR exists
    if (wdir/'OUTCAR').is_file():

        # Get content of the OUTCAR
        with open(wdir/'OUTCAR') as f: outcar_content = f.readlines()

        # Detect if VASP ended the run
        vasp_done = False
        for line in outcar_content:
            if line.find('General timing')!=-1: vasp_done = True
            if line.find('I REFUSE TO CONTINUE WITH THIS SICK JOB')!=-1: 
                status[i] = 'aborted'
                done[i] = True

        if vasp_done:
        
            # Get the global number of ionic steps
            step_tot = 0 
            step = 0
            for xf in wdir.glob("XDATCAR*"):
                with open(xf) as f:
                    for line in f:
                        if line.find('configuration=')!=-1: 
                            step = int(line.split()[-1])
                step_tot += step
            nsteps[i] = step_tot
        
            # Get the projectile distance from the surface
            with open(wdir/'XDATCAR') as f: xdatcar_content = f.readlines()
            # Get vertical size of the cell
            scale = float(xdatcar_content[1].split()[-1])
            zcell = float(xdatcar_content[4].split()[-1])*scale
            # Get total number of atoms
            nions = sum([int(s) for s in xdatcar_content[6].split()])
            for iline,line in enumerate(xdatcar_content):
                if line.find('Direct')!=-1:
                    z[i] = (float(xdatcar_content[iline+nions].split()[2]) \
                       -float(xdatcar_content[iline+nions-1].split()[2]))*zcell

            # Check if scattered
            if z[i] > z_stop:
                status[i] = 'scattered'
                done[i] = True

            # Check if the maximum number of ionic steps is reached
            elif nsteps[i] >= nsteps_max: 
                status[i] = 'stuck'
                done[i] = True

            # Check if the cluster time elapsed
            else: 
                status[i] = 'restarting_finishing_SCC'

# work until all trajectories are done or aborted
while not all(done):
    
    # Update the list of trajectories to oversee
    traj_list = []
    for d in p.iterdir():
        if d.is_dir():
            try:
                int(d.name)
                traj_list.append(d)
            except: pass
    

    # loop over not-finished trajs
    for i, wdir in enumerate(traj_list):
    
        # Get the global number of ionic steps
        step_tot = 0 
        step = 0
        for xf in wdir.glob("XDATCAR*"):
            with open(xf) as f:
                for line in f:
                    if line.find('configuration=')!=-1: 
                        step = int(line.split()[-1])
            step_tot += step
        nsteps[i] = step_tot

        # Set restart counter
        restart_counter[i] = len(list(wdir.glob('OUTCAR_*')))

        if not done[i]:
            outcar_fname  = wdir / 'OUTCAR'
            poscar_fname  = wdir / 'POSCAR'
            xdatcar_fname = wdir / 'XDATCAR'
            stopcar_fname = wdir / 'STOPCAR'
            contcar_fname = wdir / 'CONTCAR'
            oszicar_fname = wdir / 'OSZICAR'

            vasp_done = False

            # Get the projectile distance from the surface
            try:
                with open(xdatcar_fname) as f: xdatcar_content = f.readlines()
                if len(xdatcar_content) < 8: raise Exception()
            except:
                # take the initial conf if XDATCAR not ready
                with open(poscar_fname) as f: xdatcar_content = f.readlines()

            # Get vertical size of the cell
            scale = float(xdatcar_content[1].split()[-1])
            zcell = float(xdatcar_content[4].split()[-1])*scale
            # Get total number of atoms
            nions = sum([int(s) for s in xdatcar_content[6].split()])

            for iline,line in enumerate(xdatcar_content):
                if line.find('Direct')!=-1:
                    z[i] = (float(xdatcar_content[iline+nions].split()[2]) \
                       -float(xdatcar_content[iline+nions-1].split()[2]))*zcell

            # check if started
            if outcar_fname.is_file():
                # Get content of the OUTCAR
                with open(outcar_fname) as f: outcar_content = f.readlines()

                if start_time[i]==0: 
                    # Get start time
                    t_start = outcar_content[2].split()[-2]+' '+outcar_content[2].split()[-1]
                    start_time[i] = time.mktime(time.strptime(t_start, "%Y.%m.%d %H:%M:%S"))
 
                # Detect if VASP ended the run
                for line in outcar_content:
                    if line.find('General timing')!=-1: vasp_done = True
                    if line.find('I REFUSE TO CONTINUE WITH THIS SICK JOB')!=-1: 
                        status[i] = 'aborted'
                        done[i] = True
                if done[i]: continue
                        
                if vasp_done:
                
                    if status[i].find('scattered_finishing_SCC') != -1:
                    
                        # Scattering completed
                        status[i] = 'scattered'
                        done[i] = True
                        continue
                        
                    elif status[i].find('restarting_finishing_SCC') != -1:
                    
                        # Prepare for resubmission
                        status[i] = 'restarting'
                        # Save important files from the previous run
                        outcar_fname.rename( wdir/f'{outcar_fname.name}_{restart_counter[i]}')
                        poscar_fname.rename( wdir/f'{poscar_fname.name}_{restart_counter[i]}')
                        xdatcar_fname.rename(wdir/f'{xdatcar_fname.name}_{restart_counter[i]}')
                        oszicar_fname.rename(wdir/f'{oszicar_fname.name}_{restart_counter[i]}')
                        Path(wdir/'vasp.out').rename(wdir / f'vasp.out_{restart_counter[i]}')
                        Path(wdir/'vasprun.xml').rename(wdir / f'vasprun.xml_{restart_counter[i]}')
                        # reset timing
                        start_time[i] = 0
                        # rename CONTCAR to POSCAR
                        contcar_fname.rename(poscar_fname)
                        # restart calculations
                        os.chdir(wdir)
                        # Submit the job
                        if host=='local':
                            os.system(f'sh {job_script} &')
                        else:
                            os.system(f'sbatch {job_script}')
                        os.system('pwd')
                        os.chdir(p)
                        
                    elif status[i] == 'stuck_finishing_SCC':
                        status[i] = 'stuck'
                        done[i] = True

                    # Time ran out (stuck)
                    elif status[i] == 'running':
                        status[i] = 'stuck'
                        done[i] = True

                    else:
                        print(f'Error: {status[i]} here should never happen')
                        done[i] = True
                    
                else: # VASP is running
                
                    # Set timer
                    current_time[i] = time.time()

                    if stopcar_fname.is_file():
                        if current_time[i] - stopcar_fname.stat().st_ctime > max_stopcar_life:
                            status[i] = 'broken'
                            done[i] = True
                        continue

                    else:
                
                        status[i] = 'running'
                        
                        # Check if the cluster time elapsed
                        if current_time[i]-start_time[i] > max_time: 
                            # Create STOPCAR and
                            with open(stopcar_fname,'w') as f:
                                f.write('LSTOP = .TRUE.\n')
                            status[i] = 'restarting_finishing_SCC'

                        # Check if the maximum number of ionic steps is reached
                        if nsteps[i] >= nsteps_max: 
                            # Create STOPCAR and
                            with open(stopcar_fname,'w') as f:
                                f.write('LSTOP = .TRUE.\n')
                            status[i] = 'stuck_finishing_SCC'
                            
                        # Check if scattered and
                        # Create STOPCAR if scattered criterion satisfied
                        if z[i] > z_stop:
                            with open(stopcar_fname,'w') as f: 
                                f.write('LSTOP = .TRUE.\n')
                            status[i] = 'scattered_finishing_SCC'

            else: # OUTCAR is not found
                if (restart_counter[i]>0): status[i] = 'restarting'

    # Save status info into 
    with open(runchk_fname,'w') as f:
        # Create a header
        f.write(f'{"n_traj":>8}  {"status":>30} {"done":>6} {"nsteps":>6} {"rc":>6} {"start_time":>25} {"z":>8}\n')
        # Sort trajectories
        i_s = sorted(range(len(nsteps)), key=nsteps.__getitem__)
        # Write the info
        for i in i_s:
            cur_date = 0
            if start_time[i]>0: 
                cur_date = time.strftime("%d %b %Y %H:%M:%S",time.localtime(start_time[i]))
            f.write(f'{traj_list[i].name:>8}  {status[i]:>30} {done[i]:>6}' + 
                    f' {nsteps[i]:>6} {restart_counter[i]:>6}'+ 
                    f' {cur_date:>25} {z[i]:>8.2f}\n')

    # Sleep
    time.sleep(chk_time)

print('runchk says goodbye.')  
