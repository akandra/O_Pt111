#%%
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model
from scipy import special
import matplotlib.animation as animation
import functools

# Define the fitting function
def erf(x,a,x0,w):
    return a*(special.erf((x-x0)/w) + 1)/2

def animate(frame):
        
    plt.cla()

    bins = np.linspace(-4, 5, 50)

    axis.set_xlim(-3,5)
    axis.set_ylim(top=40) 
    axis.set_xlabel(r'$z_C$ / $\AA$')
    axis.set_ylabel('# per bin')
    axis.text(0.025,0.95,f'{frame} fs', transform = axis.transAxes)
    
    plt.hist(height_data[frame,:], bins=bins)


p = Path(__file__).resolve().parent
# Get data folders
if len(sys.argv)==1:
    print(f"Usage: {sys.argv[0]} [path_to_data_1] [path_to_data_2] ...")
    exit()

calc_type = sys.argv[1:]
wdirs = [p / dd for dd in calc_type]


dfs = [pd.read_fwf(d / "info.dat") for d in wdirs]


# C-atom z-position time-dependent distribution

n_frames = 100

height_data = []
for itraj in dfs[0].n_traj:

    trdir = wdirs[0] / str(itraj)
    
    try:
        height_data.append(np.loadtxt(trdir / 'height_C.dat'))
    except:
        #height_data.append(np.array([]))
        pass

df = dfs[0]
ntrajs = len(df)
height_data = -10*np.ones((n_frames,ntrajs))

for i, itraj in enumerate(df.n_traj):

    trdir = wdirs[0] / str(itraj)
    try:
        tmp = np.loadtxt(trdir / 'height_C.dat')[:,1]
        if len(tmp) < n_frames:
            height_data[:len(tmp),i] = tmp
            if df.status[i] == 'scattered':
                if df.z[i]< 4.0: height_data[len(tmp):,i] = 4.9
                else:            height_data[len(tmp):,i] = tmp[-1]
            else:
                height_data[len(tmp):,i] = tmp[-1]
        else:
            height_data[:,i] = tmp[:n_frames]
            
    except: pass

#%%
fig, axis = plt.subplots()
fig.suptitle('Position distribution for non-scattered C-atoms',
          fontweight = "bold")

anim = animation.FuncAnimation(fig, animate, frames = n_frames, interval=1,
                                    repeat=False)
 
writervideo = animation.FFMpegWriter(fps=25) 
anim.save('height_distr.mp4', writer=writervideo) 
plt.show() 

