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

p = Path(__file__).resolve().parent
calc_type = ["triplet", "unpolarized"]
wdirs = [p / dd for dd in calc_type]


dfs = [pd.read_fwf(d / "info.dat") for d in wdirs]

#%%
def animate(frame, bar_container, height_data, frame_text):
    z = []
    for height in height_data:
        try:
            z.append(height[frame,1]) 
        except: pass
    z_n, _ = np.histogram(z, HIST_BINS)
    for count, rect in zip(z_n, bar_container.patches):
        rect.set_height(count)    
    frame_text = axis.text(0.025,0.95,f'{frame} fs', transform = axis.transAxes)
    return bar_container.patches



# C-atom z-position time-dependent distribution

height_data = []
for itraj in dfs[0].n_traj:

    trdir = wdirs[0] / str(itraj)
    
    try:
        height_data.append(np.loadtxt(trdir / 'height_C.dat'))
    except:
        #height_data.append(np.array([]))
        pass

#%%
HIST_BINS = np.linspace(-5, 5, 20)

z_hist = []
for height in height_data:
    try:
        z_hist.append(height[0,1]) 
    except: pass

z_n, _ = np.histogram(z_hist, HIST_BINS)

fig, axis = plt.subplots()
fig.suptitle('Position distribution for non-scattered C-atoms',
          fontweight = "bold")

#axis.set_xlim(-3,5)
axis.set_ylim(top=30) 
axis.set_xlabel(r'$z_C$ / $\AA$')
axis.set_ylabel('# per bin')


_, _, bar_container  = axis.hist(z_hist, HIST_BINS, lw=1,
                              ec="yellow", fc="green", alpha=0.5)

frame_text = axis.text(0.025,0.95,'0 fs', transform = axis.transAxes)

anim = functools.partial(animate, bar_container=bar_container,
                                  height_data=height_data, 
                                  frame_text=frame_text)
ani = animation.FuncAnimation(fig=fig, func=anim, frames=1000, repeat=False, blit=True)
 
#plt.savefig("z_final.png") 
plt.show() 

