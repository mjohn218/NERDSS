#!/usr/bin/env python

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import sys

mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
font = {'family' : 'sans-serif',
        'size' : 22}
mpl.rc('font', **font)

df = [] 
for i in range(1, int(sys.argv[1])):
    df.append(pd.read_csv("run" + str(i) + "/species.dat"))

# get average of the specie
dfCat = pd.concat(df).groupby(level=0)
avg = dfCat.mean()
sem = dfCat.sem()
comb = avg.join(sem, lsuffix="_mn", rsuffix="_sem")

# Plot it
fig = plt.gca()

comb['Itr_mn'] = comb['Itr_mn'].apply(lambda x: x * 1E-7)
plt.fill_between(comb['Itr_mn'], comb['leg_mn'] + comb['leg_sem'], comb['leg_mn'] - comb['leg_sem'], color='k', alpha=0.25, label="_")
comb.plot(x='Itr_mn', y="leg_mn", color='k', lw=1, ax=fig, label=r"K\(_{\text{D}}\) = 0.2 \(\mu\)M")

# df = [] 
# for i in range(1, int(sys.argv[1])):
#     df.append(pd.read_csv("kd1/batch/run" + str(i) + "/species.dat"))
#
# # get average of the specie
# dfCat = pd.concat(df).groupby(level=0)
# avg = dfCat.mean()
# sem = dfCat.sem()
# comb = avg.join(sem, lsuffix="_mn", rsuffix="_sem")
#
# comb['Itr_mn'] = comb['Itr_mn'].apply(lambda x: x * 1E-7)
# plt.fill_between(comb['Itr_mn'], comb['leg_mn'] + comb['leg_sem'], comb['leg_mn'] - comb['leg_sem'], color='r', alpha=0.25, label="_")
# comb.plot(x='Itr_mn', y="leg_mn", color='r', lw=1, ax=fig, label=r"K\(_{\text{D}}\) = 1 \(\mu\)M")
#
# df = [] 
# for i in range(1, int(sys.argv[1])):
#     df.append(pd.read_csv("kd100/batch/run" + str(i) + "/species.dat"))
#
# # get average of the specie
# dfCat = pd.concat(df).groupby(level=0)
# avg = dfCat.mean()
# sem = dfCat.sem()
# comb = avg.join(sem, lsuffix="_mn", rsuffix="_sem")
#
# comb['Itr_mn'] = comb['Itr_mn'].apply(lambda x: x * 1E-7)
# plt.fill_between(comb['Itr_mn'], comb['leg_mn'] + comb['leg_sem'], comb['leg_mn'] - comb['leg_sem'], color='b', alpha=0.25, label="_")
# comb.plot(x='Itr_mn', y="leg_mn", color='b', lw=1, ax=fig, label=r"K\(_{\text{D}}\) = 100 \(\mu\)M")

fig.ticklabel_format(axis='x', style='sci', scilimits=(0,2))
# fig.set_title(r"Number of bound leg pairs, irreversible hexagon formation")
fig.set_xlabel(r"Time(s)")
fig.set_ylabel(r"Bound Legs")
fig.set_xscale('log')
fig.get_legend().remove()

plt.savefig("species.png", dpi=1000, bbox_inches='tight')
plt.close()
