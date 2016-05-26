# Adrian Soto
# May 26, 2016
# Stony Brook University
# 
# Example to calculate local order parameters of water
# from the atomic coordinates using the Waterpack package
# 
# Up until know it can only calculate dimer coordinates
# but the code structure allows any cluster size
#
#
import math
import numpy as np
import h2o



######################################################################
##############       PROGRAM    STARTS    HERE        ################
######################################################################



# Input data file
ANIfile="water_100steps.ANI"

#number of atoms per molecule
napm=3 


# Read file
nat=[]
comments=[]
txyz=[] # list of lists containing [atomtype, x,y,z]
nat, comments, txyz = h2o.readANI(ANIfile)


nsnap=len(nat)


# Lattice vectors and lattice matrix


# Simulation cell vectors and transformation matrix
ct="cartesian"
l=14.7222 # (Ang)
a1=[l,0,0]
a2=[0,l,0]
a3=[0,0,l]
A=np.vstack((a1,a2,a3)).T # Transpose is crucial!!


# Calculation input
LSIcutoff=4.50 # (Ang)


# Header for output file
header="      q              Sk             LSI(Ang^2)"
#print header




# Loop over all shapshots in file
for isnap in range(0,nsnap):
    
    snap=h2o.snapshot(nat[isnap], A, txyz[isnap], coordtype=ct)
    allmol=h2o.cluster(nat[isnap]/napm, A, snap.tx)
    allmol.H2Oindx
    allmol.FindH2Os(snap.A, snap.Ainv, cutoff=LSIcutoff)


#    allmol.printsnap( 'L1: '+ str(A.T[0][:])+' L2: '+ str(A.T[1][:]) +' L3: '+ str(A.T[2][:]) )        
    for imol in range(0, allmol.nmol):

        q, Sk=allmol.getTOPs(snap.A, snap.Ainv, imol)                  
        LSI=allmol.LSI(snap.A, snap.Ainv, imol, rcutoff=3.7) 

        print q, Sk, LSI

