#!/usr/bin/env python3
# collect trajectory info without referring to submit.log 
# %%
from pathlib import Path
import tarfile
import sys,os
import numpy as np
import time
from ase.io import read as aseread
from ase.io import write as asewrite
from ase.atom import Atom

feV = 103.6382 # amu*(A/fs)^2 to eV conversion
z_buffer = 1.5 # scattering buffer (Angstrom)
ntrap = 300    # trapping criterion in steps

# Kinetic energy in eV for a particle of mass m amu moving with velocity v A/fs
def ke(m, v):
    return 0.5*m*feV*(v @ v)

d = Path().absolute()
print(d)

trajs = []
for itraj in range(1,1001):
    wdir = d / str(itraj)
    try:    os.chdir(wdir)
    except: pass
    else:   trajs.append(itraj)
print(f' {len(trajs):>8} trajectories detected.')
if len(trajs)==0: exit()

with open(d / 'info.dat','w') as output:
    # Create a header
    output.write(f'{"n_traj":8}  {"status":12} {"n_steps":8} ' +
                f'{"z":12} {"KE":12} {"theta":12} '+
                f'{"PE_0":12} {"PE":12}  {"E_excess":12} ' +
                f'{"dEtot":12} {"magmom_0":12} {"magmom":12}\n')
    
    for itraj in trajs:
        wdir = d / str(itraj)
        try:
            os.chdir(wdir)
        except:
            try:
                with tarfile.open(d / (str(itraj) + ".tar.gz"), "r:gz") as tar: 
                    tar.extractall()
            except:
                #print(f'{wdir.name} trajectory not found')
                continue

        incar_fname   = wdir / 'INCAR'
        poscar_fname  = wdir / 'POSCAR'
        outcar_fname  = wdir / 'OUTCAR'
        xdatcar_fname = wdir / 'XDATCAR'
        oszicar_fname = wdir / 'OSZICAR'
        
        tr_status = '-'
        nsteps = 0
        z0 = 0
        zf = 0
        E0i= 0
        E0 = 0
        Einc= 0
        EK = 0
        mm0 = 0
        mm  = 0
        E0s = []
        EKs = []
        mms = []
        pos  = []
        dt = 0
        # Get time step
        with open(incar_fname) as f:
            for line in f:
                if line.find('POTIM') != -1: 
                    dt = float( line.split('#')[0].split('=')[-1] )
        if dt == 0: 
            print('timestep not defined')
            exit()

        # Get the initial positions from the 1st POSCAR
        if len(list(wdir.glob('POSCAR_*'))) > 0:
            poscar_fname = poscar_fname.parent / 'POSCAR_0'
        conf0 = aseread(poscar_fname)
        
        # Number and masses of ions
        nions = len(conf0.get_masses())
        masses = conf0.get_masses()

        # Get the initial velocities from the 1st POSCAR
        with open(poscar_fname) as f: poscar_content = f.readlines()
        for (i,line) in enumerate(poscar_content):
            if line.strip() == '': 
                idx_vel = i + 1
                break
        vel0 = [np.array(l.split(),dtype=float) for l in poscar_content[idx_vel:idx_vel+nions]]

        # Surface level and initial height and KE
        z_surf = max(conf0.positions[:-2,2])
        height0 = conf0.positions[-1,2] - z_surf
        Einc = ke(masses[-1],vel0[-1])
        
        # Get list of available XDATCARS etc. 
        xdatcar_names = sorted(list(wdir.glob('XDATCAR_*')))
        if xdatcar_fname.is_file(): xdatcar_names.append(xdatcar_fname)
        oszicar_names = sorted(list(wdir.glob('OSZICAR_*')))
        if oszicar_fname.is_file(): oszicar_names.append(oszicar_fname)

        if len(xdatcar_names) == 0:
            tr_status = 'not_started'
                
        else:
            # Get positions of the last ion
            confs = []
            for xdatcar in xdatcar_names:
                confs.extend(aseread(xdatcar,index=':'))
            [pos.append(c.positions[-1]) for c in confs]

            # Create a single XDATCAR.all file
            #asewrite(wdir / 'XDATCAR.all')
            os.chdir(wdir)
            cline = f"cat {' '.join([x.name for x in xdatcar_names])} > XDATCAR.all"
            os.system(cline)
            os.chdir(d)

            # Get pot. and kin. energies, as well as mag. moment if available
            for oszicar in oszicar_names:
                with open(oszicar) as f:
                    for line in f:
                        if line.find('E0=') != -1:
                            E0s.append(float(line.split()[8]))
                            EKs.append(float(line.split()[10]))
                            if line.find('mag=') != -1:
                                mms.append(float(line.split()[16]))
                            else:
                                mms.append(float(0))
                                
            nsteps = len(pos)
            # Initial and final heights of projectile, PE and KE of the system
            z0 = pos[ 0][2] - z_surf
            zf = pos[-1][2] - z_surf
            Epot0= E0s[ 0]
            Epot = E0s[-1]
            EKtot0 = EKs[0]
            EKtot = EKs[-1]
            # Final velocity, scattering angle and KE of the projectile
            vf = (pos[-1] - pos[-2])/dt
            costheta = vf[2]/np.linalg.norm(vf)
            theta = np.arccos(costheta)*180/np.pi
            EKf = ke(masses[-1],vf)
            # Scattering energy excess in z direction
            Eexcess = EKf*costheta*costheta + Epot - Epot0 
            # Magnetic moments
            mm0 = mms[0]
            mm  = mms[-1]

            # Define the trajectory status
            if nsteps > ntrap: tr_status = 'trapped'

            if zf > z0: 
                tr_status = 'scattered'
                # Correct final KE
                EKf = EKf + (Epot - Epot0)
            else:
                if zf > z0 - z_buffer:
                    if (vf[2] > 0) and (Eexcess > 0):
                        tr_status = 'scattered'
                        # Correct final KE
                        EKf = EKf + (Epot - Epot0)
                                        
        output.write(f'{itraj:<8}  {tr_status:<12} {nsteps:<8} ' + 
                    f'{zf:<12.4f} ' +
                    f'{EKf:<12.4f} {theta:<12.4f} ' +
                    f'{Epot0:<12.4f} {Epot:<12.4f} {Eexcess:<12.4f} ' +
                    f'{Epot+EKtot-Epot0-EKtot0:<12.4f} ' +
                    f'{mm0:<12.4f} {mm:<12.4f}\n')
        try:
            np.savetxt(wdir / "height_C.dat",
                   np.vstack((np.arange(len(pos))*dt,
                              np.array(pos)[:,2]-z_surf)).transpose(),
                   fmt='%10.1f %10.3f')
        except: pass

