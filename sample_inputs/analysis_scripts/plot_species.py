#!/usr/bin/env python3
# usage: plot_species $SPECIESFILE $SPECIES1 ...
# where $SPECIES1 ... is any number of species from the observables file

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import sys

mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
font = {'family' : 'sans-serif',
        'size' : 16}
mpl.rc('font', **font)

df = pd.read_csv(sys.argv[1])

# Plot it
fig = plt.gca()

for specie in sys.argv[2:]:
    df.plot(kind='line', x='Itr', lw=1, y=specie, ax=fig)

fig.ticklabel_format(axis='x', style='sci', scilimits=(0,2))
fig.set_xlabel(r"Timestep")
fig.set_ylabel(r"Counts")

plt.savefig("species.png", dpi=300, bbox_inches='tight')
plt.show()
