# Adrian Soto
# November 28, 2016
# Stony Brook University
# 
# Example to calculate H-bond network structure.
# 
#
import math
import numpy as np
import h2o



######################################################################
##############       PROGRAM    STARTS    HERE        ################
######################################################################



# Input data file
#ANIfile="./water_1step.ANI"
ANIfile="./water_100steps.ANI"

#number of atoms per molecule
napm=3 


# OO cutoff for H2O identification
OOcutoff = 4.5  # Ang

# Read file
nat=[]
comments=[]
txyz=[] # list of lists containing [atomtype, x,y,z]
nat, comments, txyz = h2o.readANI(ANIfile)


nsnap=len(nat)


# Simulation cell vectors and transformation matrix
ct="cartesian"
l=14.7222 # (Ang)
a1=[l,0,0]
a2=[0,l,0]
a3=[0,0,l]
A=np.vstack((a1,a2,a3)).T # Transpose is crucial!!


# Header for output file
header=" INSERT HEADER HERE"
#print header


# Loop over all shapshots in file
for isnap in range(0,nsnap):
    
    snap=h2o.snapshot(nat[isnap], A, txyz[isnap], coordtype=ct)
    allmol=h2o.cluster(nat[isnap]/napm, A, snap.tx)
    allmol.H2Oindx
    allmol.FindH2Os(snap.A, snap.Ainv, cutoff=OOcutoff)
    allmol.getD(snap.A, snap.Ainv)

    
    # Construct adjacency matrix
    Adj=np.matrix(np.zeros((allmol.nmol, allmol.nmol), dtype=np.int))
    for imol1 in range(0, allmol.nmol-1):
        for imol2 in range(imol1+1, allmol.nmol):
        
            HB=allmol.isthisHB(imol1, imol2, 1, snap.A, snap.Ainv)
            
            if (HB == 1):
                Adj[imol1,imol2] = 1
            if (HB == -1):
                Adj[imol2,imol1] = 1

    # Calculate powers of Adj
    A1=Adj.astype(int)
    A2=np.dot(A1,A1)
    A3=np.dot(A1,A2)
    A4=np.dot(A1,A3)
    A5=np.dot(A1,A4)
    A6=np.dot(A1,A5)
    A7=np.dot(A1,A6)
    A8=np.dot(A1,A7)
    A9=np.dot(A1,A8)
    A10=np.dot(A1,A9)
    A11=np.dot(A1,A10)
    A12=np.dot(A1,A11)
    A13=np.dot(A1,A12)
    A14=np.dot(A1,A13)
    A15=np.dot(A1,A14)
    #A16=np.dot(A1,A15)
    #A17=np.dot(A1,A16)
    #A18=np.dot(A1,A17)
    #A19=np.dot(A1,A18)
    #A20=np.dot(A1,A19)


    # Get number of H-bonds and loops
    for imol in range(0,allmol.nmol):
        Nacc = np.sum(Adj[:,imol])   # Number of accepting bonds 
        Ndon = np.sum(Adj[imol,:])   # Number of donating bonds 
        N3   = A3[imol,imol]         # Number of loops of length 3
        N4   = A4[imol,imol]         # Number of loops of length 4
        N5   = A5[imol,imol]         # Number of loops of length 5
        N6   = A6[imol,imol]         # Number of loops of length 6
        N7   = A7[imol,imol]         # Number of loops of length 7
        N8   = A8[imol,imol]         # Number of loops of length 8
        N9   = A9[imol,imol]         # Number of loops of length 9
        N10  =A10[imol,imol]         # Number of loops of length 10
        N11  =A11[imol,imol]         # Number of loops of length 11
        N12  =A12[imol,imol]         # Number of loops of length 12
        N13  =A13[imol,imol]         # Number of loops of length 13
        N14  =A14[imol,imol]         # Number of loops of length 14
        N15  =A15[imol,imol]         # Number of loops of length 15
        #N16  =A16[imol,imol]         # Number of loops of length 16
        #N17  =A17[imol,imol]         # Number of loops of length 17
        #N18  =A18[imol,imol]         # Number of loops of length 18
        #N19  =A19[imol,imol]         # Number of loops of length 19
        #N20  =A20[imol,imol]         # Number of loops of length 20
        
        print imol+1, Nacc, Ndon, N3, N4, N5, N6, N7, N8, N9, N10, N11, N12, N13, N14, N15 #, N16, N17, N18, N19, N20
