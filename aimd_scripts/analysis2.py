#%%
from pathlib import Path
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from lmfit import Model
from scipy import special
import matplotlib.animation as animation
import functools

# Define the fitting function
def erf(x,a,x0,w):
    return a*(special.erf((x-x0)/w) + 1)/2

p = Path(__file__).resolve().parent

# Get data folders
if len(sys.argv)==1:
    print(f"Usage: {sys.argv[0]} [path_to_data_1] [path_to_data_2] ...")
    exit()

calc_type = sys.argv[1:]
wdirs = [p / dd for dd in calc_type]


dfs = [pd.read_fwf(d / "info.dat") for d in wdirs]



#%%
# Energy loss spectra
num_bins = 20

for idf, df in enumerate(dfs):

    n0, bins0, patches0 = plt.hist(
        df[(df['status']=='scattered')].KE, 
        label='all',
        bins=num_bins,
        range=(0,1.8))
    
    n, bins, patches = plt.hist(
        df[(df['status']=='scattered') & 
        (df['magmom']>1.5) & 
        (abs(df['dEtot'])<0.1)].KE, 
        label=r'$\mu>1.5$, $|\Delta E| < 0.1$',
        bins=num_bins,
        range=(0,1.8),
        alpha = 0.7)

    plt.legend(prop ={'size': 10})

    plt.xlabel('kinetic energy / eV')
    plt.ylabel('amplitude')
    
    plt.title(f' C-atom energy distribution: {calc_type[idf]}',
            fontweight = "bold", fontsize = 12)

    plt.savefig(wdirs[idf] / "en_scattered.png") 
    plt.show() 

#%%
# Scattered trajs 

for idf, df in enumerate(dfs):
    plt.scatter(df[df['status']=='scattered'].n_steps,
            df[df['status']=='scattered'].KE,
                label='all')

    plt.scatter(df[(df['status']=='scattered') & 
        (df['magmom']>1.5) & 
        (abs(df['dEtot'])<0.1)].n_steps,
            df[(df['status']=='scattered') & 
        (df['magmom']>1.5) & 
        (abs(df['dEtot'])<0.1)].KE,
        label=r'$\mu>1.5$, $|\Delta E| < 0.1$')


    plt.xlabel('final time / fs')
    plt.ylabel('energy / eV')
    
    plt.title(f'Scattered C-atom Energy: {calc_type[idf]}',
            fontweight = "bold", fontsize = 12)
    plt.legend(prop ={'size': 10})

    plt.savefig(wdirs[idf] / "time_scattered.png") 
    plt.show()

#%%
#Scattering probability
dt = 1 # timestep in fs

for idf, df in enumerate(dfs):

    n_tot=len(df)
    n_sc=len(df[df['status']=='scattered'])
    n_tr=len(df[df['status']=='trapped'])
    n_notst=len(df[df['status']=='not_started'])
    x = np.array(range(0,350,20))
    y = [len(df[(df['n_steps']<=n) & 
                (df['status']=='scattered')])/(n_tot-n_notst) for n in x]

    # Defining the various parameters
    erf_model = Model(erf)
    fit = erf_model.fit(y, x=x, a=1, x0=200, w=1)
    print(fit.fit_report())

    # Scattered trajs 
    plt.scatter(x,y)

    a  = fit.params['a'].value
    x0 = fit.params['x0'].value
    w  = fit.params['w'].value
    x1 = range(0,1000,10)
    #plt.plot(x1, erf(x1,a,x0,w),c='red',ls='-',lw=2)

    plt.xlabel('time / fs')
    plt.ylabel('probability')
    
    plt.title(f'Estimator for scattering probability: {calc_type[idf]}',
            fontweight = "bold", fontsize = 12)
    #plt.legend(prop ={'size': 10})

    plt.savefig(wdirs[idf] / "scattering_probability.png") 
    plt.show()

#%%
#Scattering probability combined
dt = 1 # timestep in fs

#cmap = mpl.colormaps['viridis']
#with plt.style.context("grayscale"):

fig, axis = plt.subplots()
fig.suptitle(f'Estimator for scattering probability',
            fontweight = "bold", fontsize = 12)

if len(dfs)>1:
    for idf, df in enumerate(dfs):

        n_tot=len(df)
        n_sc=len(df[df['status']=='scattered'])
        n_tr=len(df[df['status']=='trapped'])
        n_notst=len(df[df['status']=='not_started'])
        x = np.array(range(0,410,20))
        y = [len(df[(df['n_steps']<=n) & 
                    (df['status']=='scattered')])/(n_tot-n_notst) for n in x]

        # Scattered trajs 
        axis.scatter(x,y, label=calc_type[idf])

    plt.xlabel('time / fs')
    plt.ylabel('probability')
        
    plt.legend(prop ={'size': 10})

    plt.savefig("scattering_probability.png") 
    plt.show()


