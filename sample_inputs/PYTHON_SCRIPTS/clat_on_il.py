#!/usr/bin/env python3                                                                                                                                            
# usage: clat_on_mem $complex_filename StringName                                                                                                                               
# finds all complexes that include clathrin as the third element (pip, ap, clat, *)
#for each iter, store the n clathrin in each of these cluster, sorts, and prints the largest cluster.
# outfile: MaxClathCluster.txt printes out [iteration Ncluster MaxNClathrinPerCluster]

#import matplotlib as mpl                                                                                                                                         
#import matplotlib.pyplot as plt                                                                                                                                  
#import pandas as pd                                                                                                                                              
import sys
import os





filepath = sys.argv[1]
srank=sys.argv[2]

fname="ArrayClathrinMem%s.txt" % (srank)
fname3="ClathrinOnMem%s.txt" % (srank)
fname2="MaxClathCluster%s.txt" % (srank)

outfile=open(fname,"w")
outfile2=open(fname2,"w")
outfile3=open(fname3,"w")

print("File path {}: ".format(filepath))
if not os.path.isfile(filepath):
    print("File path {} does not exist. Exiting...".format(filepath))
    sys.exit()

#inum=[0]*10000
#cnum=[0]*100000
#comnum=[0]*100000
order=[]

numbers=[1, 3, 4, 2]
numbers.sort()
print(numbers)
justnum='0'
with open(filepath) as fp:
    icnt = 0
    ccnt = 0 
    newit=0
    sum = 0
    for line in fp:
#            print("line {} contents {}".format(icnt, line))                                                                                                          
        
        if(line.find('iter') != -1):
            #sort the order list which is of length newit
            ni=newit
            print("Print from last iter: {} {}".format(justnum[0], newit))
            
            if(ni>0):
                #sort the elements of order
                order.sort()
                
                outfile2.write("{} {} {}\n".format(justnum[0], newit, order[newit-1]))#previous iteration, and last element sorted
                #print(order)
                #for i in range(0, ni):
                #    outfile2.write("e: {} ".format(order[i]))
                #outfile2.write("\n")
            else:
                outfile2.write("{} {} {}\n".format(justnum[0], newit, 0))#previous iteration, and last element sorted
                               
            outfile.write("{}\n".format(justnum[0]))
            outfile3.write("{} {}\n".format(justnum[0], sum))
                
            arr1 = line.split(' ')
            justnum = arr1[1].split('\n')
            newit=0
            del order[:]

#            print(len(justnum))
#            print(justnum[0])
            #inum[icnt]=justnum[0]
            #cnum[icnt]=0
            sum = 0            
            #print("iter")
            #print(justnum[0])                                                                                                                                    
            icnt +=1
        if(line.find('clat:') !=-1):
           arr=line.split('\t')
           #cnum[icnt-1]=arr[0]
           
#           cname=arr[1].split('\n')
           values=arr[1].split(' ')
           #print("read in")
           #print(arr[1])
           lt=len(values)
           #print(len(values))
           #print(values[0])
           #print(values[lt-2])
           if(len(values)>4):   #these are clat: n IL: m
                if(values[2] == 'IL:'):
                    #print(values[5])
                    nclatarr=values[1].split('.')
                    nclat=nclatarr[0]
                    #comnum[ccnt]=nclat
                    outfile.write("{} {}\n".format(arr[0],nclat))
                    sum+= int(nclat)*int(arr[0])
                    ccnt+=1
                    order.append(int(nclat))#[newit]=nclat
                    newit+=1
           else:
               #this is monomeric clathrin, above has dimers and mem bound
               nclatarr=values[1].split('.')
               nclat=nclatarr[0]
                   #comnum[ccnt]=nclat
               #outfile.write("{} {}\n".format(arr[0],nclat))
               ccnt+=1
               order.append(int(nclat))#[newit]=nclat
               newit+=1
                    


#last line
order.sort()
if(newit>0):
    outfile2.write("{} {} {}\n".format(justnum[0], newit, order[newit-1]))#previous iteration, and last element sorted
else:
    outfile2.write("{} {} {}\n".format(justnum[0], newit, 0))#previous iteration, and last element sorted
      
outfile3.write("{} {}\n".format(justnum[0], sum))            
#print(" 2 columns!")                                                                                                                                             
#for i, iters in enumerate(inum):
#    value=cnum[i]
#    step=inum[i]
#    outfile.write("{} {}\n".format(step,value))

