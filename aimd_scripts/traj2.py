#%%
from pathlib import Path
import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
import ase.build
from ase.visualize import view
from ase.io import Trajectory
import copy

feV = 103.6382 # amu*(A/fs)^2 to eV conversion

# Kinetic energy in eV for a particle of mass m amu moving with velocity v A/fs
def ke(m, v):
    return 0.5*m*feV*(v @ v)

def read_file(fname):
    try:
        with open(fname) as f: data = f.readlines()
    except:
        print(f' {fname} not found')
        exit()
    else:
        return data

def get_conf_info(traj_dir):
    # will always return cartesian coordinates!

    xdatcar_names = sorted(list(traj_dir.glob('XDATCAR_*')))
    if (traj_dir/'XDATCAR').is_file(): 
        xdatcar_names.append(traj_dir/'XDATCAR')
    
    pos = []
    for xf in xdatcar_names:
        data    = read_file(xf)
        scale   = float(data[1].split()[-1])
        a1      = np.array(data[2].split(),dtype=float)*scale
        a2      = np.array(data[3].split(),dtype=float)*scale
        a3      = np.array(data[4].split(),dtype=float)*scale
        cell    = np.array([a1,a2,a3])
        # Get species and number of atoms
        species = data[5].split()
        nions   = [int(s) for s in data[6].split()]

        r = np.zeros((sum(nions),3))
        for il,line in enumerate(data):
            if "Direct configuration=" in line:
                for j in range(sum(nions)):
                    r[j,:] = cell @ np.array(data[il+j+1].split()[:3],dtype=float)
                pos.append(copy.copy(r))

    return pos

def get_energy_info(traj_dir):
    # will always return cartesian coordinates!

    oszicar_names = sorted(list(traj_dir.glob('OSZICAR_*')))
    if (traj_dir/'OSZICAR').is_file(): 
        oszicar_names.append(traj_dir/'OSZICAR')
    
    pe = []
    mm = []
    ket= []
    for xf in oszicar_names:
        with open(xf) as f:
            for line in f:
                if 'E0=' in line:
                    pe.append(float(line.split()[8]))
                    ket.append(float(line.split()[10]))
                    if 'mag=' in line:
                        mm.append(float(line.split()[16]))
                    else:
                        mm.append(float(0))
    return pe, mm, ket

p = Path(__file__).resolve().parent

#df = pd.read_fwf("info.dat")

# Get trajectory number
if len(sys.argv)==1:
    print(f"Usage:{sys.argv[0]} [path_to_traj_1] [path_to_traj_2] ...")
    exit()

traj_dirs = []
for arg in sys.argv[1:]:
    traj_dirs.append(p / arg)

#traj_dir = p / '41'
#results_dir = p.mkdir(exist_ok=True)

# Consider to get the values from INCARs
Einc = 1.46*np.ones(len(sys.argv)-1)
dt   = 1.00*np.ones(len(sys.argv)-1)
m = 12

pes = []
mms = []
kets = []
rs = []
kes = []
z_surf = []

for it,traj_dir in enumerate(traj_dirs):
    confs = get_conf_info(traj_dir)
    z_surf.append(max(confs[0][:-1,2]))
    
    pe, mm, ket = get_energy_info(traj_dir)
    pes.append(pe)
    mms.append(mm)
    kets.append(ket)

    r = [ p[-1,:] for p in confs]
    vel = [ (r[i+1]-r[i])/dt[it] for i in range(len(r)-1)]
    en = [ ke(m,v) for v in vel ]
    en.insert(0,Einc[it])
    kes.append(en)
    rs.append(r)
    

fig, axes = plt.subplots(4,sharex=True,figsize=(8,10))

for i,r in enumerate(rs):
    axes[0].plot(range(len(r)),[pos[2] - z_surf[i] for pos in r],
                 label=sys.argv[i+1])

axes[0].set_xlabel('time / fs')
axes[0].set_ylabel('height / A')
axes[0].legend()

for i,pe in enumerate(pes):
    axes[1].plot(range(len(pe)),[e-pes[0][0] for e in pe],
    label=f'pe {sys.argv[i+1]}')
for i,ken in enumerate(kes):
    axes[1].plot(range(len(ken)),[k-ken[0] for k in ken],ls="dotted",
    label=f'ke {sys.argv[i+1]}')

axes[1].set_ylim(-4,6)
axes[1].set_xlabel('time / fs')
axes[1].set_ylabel(r'$\Delta E$ / eV')
axes[1].legend()
if len(pes)>1: 
    axes[1].text(0.01,0.95,rf'$\Delta E=${(pes[0][0] - pes[1][0]):4.2f}eV',
     transform = axes[1].transAxes)

for i,mm in enumerate(mms):
    axes[2].plot(range(len(mm)),[np.abs(m) for m in mm])
axes[2].set_xlabel('time / fs')
axes[2].set_ylabel(r'$|N_\uparrow$ - $N_\downarrow|$')

for i,pe in enumerate(pes):
    axes[3].plot(range(len(pe)),[kets[i][t] + pe[t] - kets[i][0] -pe[0] for t in range(len(pe))],
    label=f'{sys.argv[i+1]}')

#axes[3].set_ylim(-4,6)
axes[3].set_xlabel('time / fs')
axes[3].set_ylabel(r'$\Delta E_t$ / eV')
axes[3].legend()

for ax in axes:
    ax.label_outer()
    

fig.tight_layout()
plt.subplots_adjust(hspace=0.0)
plt.savefig(f"traj_147.png") 
plt.show()

exit()
temp = ase.io.read(traj_dir / 'vasprun.xml')
view(temp)

#%%
# Energy loss spectra
num_bins = 5
  
n, bins, patches = plt.hist(
    df[(df['status']=='scattered') & 
       (df['magmom']>1.5) & 
       (abs(df['dEtot'])<0.1)].KE, 
    bins=num_bins,
    range=(0,1.8))
  
plt.rcParams.update({'font.size': 14})

plt.xlabel('kinetic energy / eV')
plt.ylabel('amplitude')
 
plt.title('Final C-atom energy',
          fontweight = "bold")

plt.savefig("en_scattered.png") 
plt.show() 

#%%
# Scattered trajs 
plt.scatter(df[df['status']=='scattered'].n_steps,
         df[df['status']=='scattered'].KE)

plt.xlabel('time / fs')
plt.ylabel('energy / eV')
 
plt.title('Energy of scattered C-atoms vs. final time',
          fontweight = "bold")

plt.savefig("time_scattered.png") 

#%%
#Sticking probability
n_tot=len(df)
n_sc=len(df[df['status']=='scattered'])
n_tr=len(df[df['status']=='trapped'])
n_notst=len(df[df['status']=='not_started'])
st_p_upper = 1 - n_sc/(n_tot-n_notst)
st_p_lower = n_tr/(n_tot-n_notst)
[len(df[df['n_steps']>n])/(n_tot-n_notst) for n in range(300,1100,50)]

#%%
# Final C-atom z-position
num_bins = 20
  
n, bins, patches = plt.hist(
    df[df['status']!='scattered'].z, 
    bins=num_bins,
    range=(-2,3))
  
plt.xlabel(r'$z_C$ / $\AA$')
plt.ylabel('# per bin')
 
plt.title('Position distribution for non-scattered C-atoms',
          fontweight = "bold")

plt.savefig("z_final.png") 
plt.show() 
