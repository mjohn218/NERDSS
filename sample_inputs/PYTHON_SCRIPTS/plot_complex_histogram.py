#!/usr/bin/env python2
# Created 06/06/19, MJV
# Plots the complex component histogram of the final frame
# ONLY WORKS FOR HOMOGENOUS COMPLEXES FROM HOMOGENOUS SIMULATIONS

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)
import pandas as pd
import argparse
import numpy as np
import seaborn as sns

import matplotlib as mpl
font = {'family' : 'sans-serif',
        'size' : 20 }
mpl.rc('font', **font)

parser = argparse.ArgumentParser(description="Plots the histogram of complex sizes.")
reqdArgs = parser.add_argument_group('required arguments')
optArgs = parser.add_argument_group('optional arguments')

# set up arguments
reqdArgs.add_argument('-f', action='append', dest='file', type=str, help='Complex histogram file, output from NERDSS as a file named ComplexHistogram_Np$NUMPROS.dat')
optArgs.add_argument('-x', '--xlabel', dest='xlabel', type=str, help='X-axis label', default='Species per complex')
optArgs.add_argument('-y', '--ylabel', dest='ylabel', type=str, help='Y-axis label', default='Counts')
optArgs.add_argument('-t', '--title', dest='title', type=str, help='Title of the plot', default=None)
args = parser.parse_args()

sns.set(style="ticks", font_scale=1.5, context='paper')

if len(args.file) == 1:
    rawData = [line.split() for line in open(args.file[0])]
    data = []
    header = "iter:"
    for elem in reversed(rawData):
        if (elem[0] != header):
            data.append([float(elem[0]), elem[1], float(elem[2])])
        else:
            break

    data.sort(key=lambda x: x[2]) # sort by the number of species in each complex
    data = zip(*data) # transpose the array

    fig, ax1 = plt.subplots()
    ax1.bar(data[2], data[0])
    # ax1.set_xticks(np.arange(1, data[2] + 1)) # make sure the x axis ticks are only integers

    data2 = []
    for i in range(0, len(data[0])):
        data2.append(data[0][i] * data[2][i])
    ax2 = plt.axes([0,0,1,1])
    ip = InsetPosition(ax1, [0.45,0.45,0.5,0.5])
    ax2.set_axes_locator(ip)
    ax2.bar(data[2], data2)
    # ax2.set_xticks(np.arange(1, data[2] + 1)) # make sure the x axis ticks are only integers

    ax1.set_xlabel(args.xlabel)
    ax1.set_ylabel(args.ylabel)
    ax2.set_ylabel("Counts * complex size")
    if args.title:
        plt.title(args.title)
    plt.savefig('complex_histogram.png', dpi=300, bbox_inches='tight')
else:
    rawData = []
    for fileName in args.file:
        rawData.append([line.split() for line in open(fileName)])
    data = []
    maxIndex = 0
    for frame in rawData:
        tmpData = []
        for elem in reversed(frame):
            if elem[0] != "iter:":
                tmpData.append([float(elem[2]), float(elem[0])])
            else:
                break
        tmpData.sort(key=lambda x: x[0])
        tmpData = map(list, zip(*tmpData))
        if max(tmpData[0]) > maxIndex:
            maxIndex = max(tmpData[0])
        data.append(pd.DataFrame(tmpData[1], index=tmpData[0]))

    for i in range(0, len(data)):
        data[i] = data[i].reindex(np.arange(1, int(maxIndex) + 1))
        data[i] = data[i].fillna(0)

    # Create a data frame containing the average values and SEM of each
    # complex size
    dfCat = pd.concat(data).groupby(level=0)
    avg = dfCat.mean()
    sem = dfCat.sem()
    comb = avg.join(sem, lsuffix="-mn", rsuffix="-sem").fillna(0)

    indices = comb.index.tolist()
    if (max(indices) % 2 != 0):
        newIndices = indices + [max(indices) + 1]
        tmpdf = pd.DataFrame([[0,0]], columns=['0-mn','0-sem'])
        # comb = comb.append(tmpdf, ignore_index=True)
        comb = comb.append(tmpdf)
        comb = comb.reindex(newIndices)
    comb = comb.fillna(0)
    indices = comb.index.tolist()
    # create the figure
    fig = plt.figure(figsize=(6,5))
    ax1 = fig.add_subplot(111)
    ax1.bar(indices, comb['0-mn'], yerr=comb['0-sem'], capsize=4)
    # ax1.set_xticks(np.arange(1, comb.index.tolist()[-1] + 1)) # make sure the x axis ticks are only integers
    ax1.xaxis.set_major_locator(mpl.ticker.FixedLocator(np.arange(1, max(indices), max(indices)/4)))
    ax1.xaxis.set_minor_locator(mpl.ticker.FixedLocator(np.arange(1, max(indices))))

    comb2 = comb.copy()
    comb2['0-mn'] = comb2['0-mn'].multiply(indices);
    comb2['0-sem'] = comb2['0-sem'].multiply(indices);
    ax2 = plt.axes([0,0,1,1])
    ip = InsetPosition(ax1, [0.35,0.35,0.6,0.6])
    ax2.set_axes_locator(ip)
    ax2.bar(indices, comb2['0-mn'], yerr=comb2['0-sem'], error_kw=dict(lw=1, capsize=2, capthick=1))

    ax2.xaxis.set_major_locator(mpl.ticker.FixedLocator(np.arange(1, max(indices), max(indices)/4)))
    ax2.xaxis.set_minor_locator(mpl.ticker.FixedLocator(np.arange(1, max(indices))))
    ax2.yaxis.set_major_locator(mpl.ticker.FixedLocator(ax1.get_yticks()))
    emptyStringLabels = ['']*(len(indices)/2)
    ax2.set_xticklabels(emptyStringLabels)
    emptyStringLabels = ['']*(len(ax1.get_yticks())/2)
    ax2.set_yticklabels(emptyStringLabels)

    ax1.set_xlabel(args.xlabel)
    ax1.set_ylabel(args.ylabel)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    if args.title:
        plt.title(args.title)
    plt.savefig('complex_histogram.png', dpi=300, bbox_inches='tight')
    plt.close()
