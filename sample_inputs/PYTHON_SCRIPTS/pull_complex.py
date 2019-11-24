#!/usr/bin/env python3                                                                                                                                            
# usage: pull_complex $complex_name                                                                                                                               
# where $complex_name must exist in ComplexHistogram.dat                                                                                                          

#import matplotlib as mpl                                                                                                                                         
#import matplotlib.pyplot as plt                                                                                                                                  
#import pandas as pd                                                                                                                                              
import sys
import os





filepath = sys.argv[1]
srank=sys.argv[2]
rank=int(srank)
fname="bm_rank%d.dat" % (rank)

outfile=open(fname,"w")

print("File path {}: ".format(filepath))
if not os.path.isfile(filepath):
    print("File path {} does not exist. Exiting...".format(filepath))
    sys.exit()

inum={}
cnum={}
with open(filepath) as fp:
    icnt = 0
    for line in fp:
#        print("line {} contents {}".format(icnt, line))                                                                                                          

        if(line.find('iter') != -1):
            arr=line.split(' ')
            justnum=arr[1].split('\n')
            inum[icnt]=justnum[0]
            cnum[icnt]=0
#            print(justnum[0])                                                                                                                                    
            icnt +=1
        if(line.find('B: 1. M: 1.') !=-1):
           arr=line.split('\t')
           cnum[icnt-1]=arr[0]
#           print(arr[0])                                                                                                                                         

#print(" 2 columns!")                                                                                                                                             
for i, iters in enumerate(inum):
    value=cnum[i]
    step=inum[i]
    outfile.write("{} {}\n".format(step,value))

