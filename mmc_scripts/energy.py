# %%
#Plot the energy for a single trajectory

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# path to the python file
p = Path.cwd()

nargs = len(sys.argv)
if nargs < 2:
    print("I need file names.")

else:
    fig, axis = plt.subplots(2,sharex=True)
    #fig.suptitle('mmc trajectory')
    axis[0].set_ylabel("energy (eV)")
    axis[1].set_ylabel("# of adsorbates")

    for enfile in [p / n for n in sys.argv[1:]]:

        with open(enfile,'r') as f:
            step,energy,nads = np.loadtxt(f,skiprows=1).transpose()

        axis[0].plot(step, energy)
        axis[1].plot(step, nads)

    #axis[0].set_xlim([0,1000000])
    #axis[0].set_ylim([max(energy)-20,max(energy)+1])
    #axis[1].set_ylim([max(nads)-250,max(nads)+1])

    for ax in axis:
        ax.set_xlabel(r'mmc step')
        ax.label_outer()
    fig.tight_layout()

    plt.show()

