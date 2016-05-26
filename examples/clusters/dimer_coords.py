# Adrian Soto
# May 26, 2016
# Stony Brook University
# 
# Example to calculate dimer coordinates using
# the Waterpack package
# 
#
import math
import numpy as np
import h2o


######################################################################
##############       PROGRAM    STARTS    HERE        ################
######################################################################



# Input data file
ANIfile="dimers.xyz"

#number of atoms per molecule
napm=3 


# Read file
nat=[]
comments=[]
txyz=[] # list of lists containing [atomtype, x,y,z]
nat, comments, txyz = h2o.readANI(ANIfile)

# Number of dimers (or snapshots)
nsnap=len(nat)


# cell vectors and matrix
ct="cartesian"
l=20.00 
a1=[l,0,0]
a2=[0,l,0]
a3=[0,0,l]
A=np.vstack((a1,a2,a3)).T # Transpose is crucial!!


# Define and print file header
header="rOO \t phi \t theta \t rOH11 \t rOH12 \t rOH21 \t rOH22 \t HOH1 \t HOH2 \t alpha \t beta \t gamma \t nu \t mu \t OdHOa"    
print header


# Loop over all shapshots in file
for icl in range(0,nsnap):
    
    # Create instance of the class snapshot
    snap=h2o.snapshot(nat[icl], A, txyz[icl], coordtype=ct)

    # Create instance of the class cluster, in this
    # case containing all molecules in the snapshot
    mycluster=h2o.cluster(nat[icl]/napm, A, snap.tx)
    mycluster.H2Oindx
    mycluster.wrap()    
    mycluster.FindH2Os(snap.A, snap.Ainv, cutoff=4.5)    


    # Find molecule roles in the dimer:
    # - centralmol labels the donor molecule (can take values 0,1) 
    # - secondarymol labels the acceptor molecule (can take values 0,1) 
    # - primaryh labels the donated H 
    centralmol, secondarymol, primaryh, secondaryh = mycluster.findHBandsort()


    # Center and orient with respect to molecule with index centralmol
    mycluster.CaO(centralmol)


    # Bring secondary molecule to positive z semispace 
    # and get 12 dimer coordinates
    mycluster.molabovexy(secondarymol)
    xdimer=mycluster.dimercoords(centralmol, secondarymol)


    # Calculate additional dimer coordinates mu, nu and d(OA, H)
    nu, mu, dODH, dOAH, ODHOA = mycluster.HBcoords(centralmol, secondarymol, primaryh)

    # Output data
    h2o.printdimerx(xdimer + [nu, mu, dODH, dOAH, ODHOA*180.0/3.141592]) # angles in deg
    ###printdimerx(xdimer + [nu, mu, ODHOA]) # angles in rad
