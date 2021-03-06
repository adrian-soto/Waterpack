# Adrian Soto
# February 28, 2016
# Stony Brook University
# 
# Program to generate relative coordinates and 
# order paramterers of water from an ANI (or XYZ) file
# 
# Up until know it can only calculate dimer coordinates
# but the code structure allows any cluster size
#
#
import math
import numpy as np



class snapshot:
    def __init__(self, Natoms, cellA, typeandcoords, coordtype='cartesian'):

        # Number of atoms in snapshot
        if (Natoms < 1):
            print "ERROR: number of atoms must be larger than 0. EXIT"
            exit()
        else:
            self.nat = Natoms        


        # Cell vectors matrix and inverse.
        # Notice that cell vectors must be arranged
        # in columns: A=[a1, a2, a3]
        if ( abs(np.linalg.det(cellA)) < 1.0E-6):
            print "ERROR: number of atoms must be larger than 0. EXIT"
            exit()
        else:
            self.A=cellA
            self.Ainv=np.linalg.inv(cellA)
        
            
        
        # Atomic type and positions. Must be a list of lists, 
        # each sublist containing [type,x,y,z] for each atom.
        # The elements are of type STRING
        if (len(typeandcoords) != Natoms):
            print "ERROR: list must have length Natoms. EXIT"
            exit()

        else:
            # Create first list/array row. Then append/stack
            self.element=[typeandcoords[0][0]]
            self.x=np.array([float(typeandcoords[0][1]),float(typeandcoords[0][2]),float(typeandcoords[0][3])])
            for iat in range(1, self.nat):
                self.element.append(typeandcoords[iat][0])
                xnew=np.array([
                        float(typeandcoords[iat][1]),float(typeandcoords[iat][2]),float(typeandcoords[iat][3])])
                self.x=np.vstack( (self.x, xnew) )


        if (coordtype=='cartesian' or coordtype=='Cartesian'):
            pass
        
        elif (coordtype=='crystal' or coordtype=='Crystal'):
            # Convert to cartesian
            self.x=np.dot(A,self.x.T).T
        else:
            print "ERROR: unrecognized coordtype. EXIT"
            exit()
            


        # Create tx. Mostly to interface easily (but not efficiently) with class cluster
        self.tx=[]
        for iat in range(0, self.nat):
            self.tx.append([self.element[iat], self.x[iat][0], self.x[iat][1], self.x[iat][2]])
        



    def printsnap(self, comment):
        #
        # Print snapshot in xyz format
        #
        
        precision=5
        
        fmt ="% "+str(precision)+"e" # will print a space for positive and a "-" for negative numbers
        #fmt="%."+str(precision)+"e" # will misalign if there are negative numbers
        #fmt="%."+str(precision)+"e" # will print + or - sign, keeping all numbers aligned


        print "\t", self.nat
        print comment
        for iat in range(0, self.nat):
            print self.element[iat], "\t", fmt % self.x[iat][0], "\t", fmt % self.x[iat][1], "\t", fmt % self.x[iat][2]

        return

                                 

    def d_PBC(self, i, j):
    #
    # Cartesian distance between 2 atoms
    # in a system with periodic boundary
    # conditions (PBC)
    #
    # x1, x2: numpy arrays containing cartesian coordinates
    # A: numpy matrix containing cell vectors
    # Ainv: inverse of A
    #
    # Steps:
    #
    # 1) calculate vector linking xi and xj, rij
    # 2) change to cell basis
    # 3) translate to mimimum image
    # 4) change to cartesian basis
    # 5) calculate euclidean distance
    # 
    
        sij=np.dot(self.Ainv, self.x[i]-self.x[j])
        sij[:]=sij[:] - np.round(sij[:])
        rij=np.dot(self.A,sij)
        return math.sqrt( np.dot(rij,rij) )


                             


class cluster:    
    #############################
    # class cluster variables
    #############################
    def __init__(self, Nmolecules, cellA, typeandcoords):
        

        # Number of molecules in cluster
        if (Nmolecules < 1):
            print "ERROR: number of molecules must be larger than 0."
        else:
            self.nmol = Nmolecules 
            self.nat = 3*self.nmol
        

        # Matrix containing lattice vectors
        if (abs(np.linalg.det(cellA)) < 1.0E-06):
            print "ERROR: lattice vectors are NOT linearly independent. EXIT"
            exit()
        else:
            self.A = cellA
            self.Ainv=np.linalg.inv(cellA)


        # Atomic type and positions. Must be a list of lists, 
        # each sublist containing [type,x,y,z] for each atom.
        # The elements are of type STRING
        if (len(typeandcoords) != 3*Nmolecules):
            print "ERROR: list alat must have length 3*number_of_molecules."
        else:
            self.tx = typeandcoords

        # All atomic coordinates.
        # The elements are a numpy.array of FLOAT
        self.x=[]
        for iat in range(0, self.nat):
            self.x.append(np.array([
                        float(self.tx[iat][1]),float(self.tx[iat][2]),float(self.tx[iat][3]) ]))             


        # Index lists for Oxygens, Hydrogens and H2Os
        #
        self.Oindx=[]
        self.Hindx=[]
        self.H2Oindx=[]

        
        # Arrays for interatomic distances. 
        # Initialize to False, since they will only 
        # be calculated if needed.
        self.rOO=False

        self.D=False


        # Some constants
        self.deg2rad=math.pi/180.0
        self.rad2deg=180.0/math.pi
        

        

    #############################
    # class cluster functions
    ############################# 

    def d(self, x1, x2):
        #
        # Cartesian distance between points x1 and x2.
        #

        return math.sqrt( (x1[0]-x2[0])**2 + (x1[1]-x2[1])**2 + (x1[2]-x2[2])**2 )



    def d_PBC(self, A, Ainv, i, j):
    #
    # Cartesian distance between 2 atoms
    # in a system with periodic boundary
    # conditions (PBC)
    #
    #  - i, j: atom indices. x[i], x[j] are 3d
    #          np.arrays of cartesian coordinates
    #  - A: numpy matrix containing cell vectors
    #  - Ainv: inverse of A
    #
    # Steps:
    #
    # 1) calculate vector linking xi and xj, rij
    # 2) change to cell basis
    # 3) translate to mimimum image
    # 4) change to cartesian basis
    # 5) calculate euclidean distance
    # 
    
        
        sij=np.dot(Ainv, self.x[i]-self.x[j] )
        sij[:]=sij[:] - np.round(sij[:])
        rij=np.dot(A, sij)
        return np.sqrt( np.dot(rij,rij) )




    def T(self, R):
        #
        # Translate by vector R.
        # R must be a numpy array
        #
        if (R.shape[0] != 3):
            print "ERROR: translation array R has wrong properties"
        else:
            for iat in range(0, self.nat):
                self.x[iat] = self.x[iat] + R
        
        return
 



    def checkifwrapped(self):
        #
        # Check if coordinates exceed +-0.5*A[i]  
        #
        # Return 1 if wrapped, 0 if not wrapped
        #
        
        status=1
        
        min=-0.5
        max= 0.5
 
        for iat in range(0, self.nat):

            # Position vector in basis where the simulation box 
            # is a cube of unit length
            sx=np.dot(self.Ainv, self.x[iat])
                        
            if (sx[0] > max or sx[0] < min):
                status=0
                break
            
            if (sx[1] > max or sx[1] < min):
                status=0
                break
            
            if (sx[2] > max or sx[2] < min):
                status=0
                break

        return status




    def wrap(self, verbose=False):
        #
        # Translate coordinates so that they are centered around (0,0,0)
        #
        
        # Scan all coordinates and see if any exceeds half alat. 

        # Minimum and maximum values of coordinates in the cell
        # basis, where the simulation cell is a unit length cube.
        min=-0.5
        max=0.5
                        
        zero=np.array([0.0, 0.0, 0.0])
        
        status=self.checkifwrapped()

        if (status != 1 and verbose==True):
            print " " 
            print "For \vec r=\sum_i alpha_i \vec a_i, wrapping coordinates to a box with"
            print min, "< alpha1 <=", max
            print min, "< alpha2 <=", max
            print min, "< alpha3 <=", max
            print " " 


            count=1 # DEBUG

            # Wrap coordinates. Loop over x,y,z coordinates and
            # translate if coordinates exceed box.
            while status != 1:

                R=zero
                
                # Find translation vector
                
                # x a1
                for iat in range(0, self.nat):

                    # Position vector in basis where the simulation box 
                    # is a cube of unit length
                    sx=np.dot(self.Ainv, self.x[iat])
                    
                    
                    if (sx[0] > max):
                        R=R-self.A[:,0]
                        break
                    elif (sx[0] <= min):
                        R=R+self.A[:,0]
                        break
                    
                # y a2
                for iat in range(0, self.nat):
                    if (sx[1] > max):
                        R=R-self.A[:,1]
                        break
                    elif (self.x[iat][1] <= min):
                        R=R+self.A[:,1]
                        break
                    
                # z a3
                for iat in range(0, self.nat):
                    if (sx[2] > max):
                        R=R-self.A[:,2]
                        break
                    elif (sx[2] <= min):
                        R=R+self.A[:,2]
                        break

            
                if(not np.array_equal(R,zero)):
                    print " The current coordinates are"
                    for iat in range(0, self.nat):
                        print self.x[iat]
                
                    # Translate
                    print "Translating by R=" , R
                    self.T(R)
                
                    print " The new coordinates are"
                    for iat in range(0, self.nat):
                        print self.x[iat]
                    

                status=self.checkifwrapped()
            
        return 



    def FindAtoms(self):
        # Loop over all atoms. Everytime an O is found,
        # scan all Hs and find the two closest.
        
        # Initialize distances to some huge number,
        # d1<=d2
        d1=1.0E+30 
        d2=1.1E+30
        

        # Atom type labels
        Olabel='O'
        Hlabel='H'

        # Indices containing Oxygens and Hydrogens
        Oind=[]
        Hind=[]
        
        err=0
        # Find indices for Oxygens and Hydrogens
        for iat in range(0, self.nat):
            if (self.tx[iat][0] == Olabel):
                Oind.append(iat)
            elif(self.tx[iat][0] == Hlabel):
                Hind.append(iat)
            else:
                print "ERROR: Atom label with index ", iat, " is not O nor H."
                err=1
                break
                
        nO=len(Oind)
        nH=len(Hind)
        #print "DEBUG-- ", nO, " Oxygens found with indices ", Oind
        #print "DEBUG-- ", nH, " Hydrogens found with indices ", Hind

        
        self.Oindx=Oind
        self.Hindx=Hind

        # return 0 if success, 1 if error
        return err



    def FindH2Os(self, A, Ainv, cutoff=3.0):
        #
        # Make index list for atoms forming H2O molecules. 
        # The criterion to make a molecule is based solely
        # in the O-H distances. 
        #
        #
        # cutoff (same units as coordinates): leaves outside 
        # the possible O-H bonds any pair with distance beyond
        # the cutoff value. For computational efficiency when
        # large cells with many atoms
        #
        err=self.FindAtoms()

        if (len(self.Oindx)<1 or len(self.Hindx)<2):
            print "ERROR: index lists for Os and Hs don't have correct length. Possibly FindAtoms() has not been called."

        elif (err==0):

            for iO in self.Oindx:            
                
                # Find all possible O-H distances for the current O
                indices=[] 
                distances=[]
                for iH in self.Hindx:

                    daux=self.d_PBC(A, Ainv, iO , iH)                    
                    

                    if (daux < cutoff):
                        indices.append(iH)
                        distances.append(daux)
                

                nH=2 # two Hydrogens per molecule
                hydrogens=[]
                # Form molecule by assigning H's with smallest O-H distance
                #for iH in range(0, alen(self.Hindx)):
                for iH in range(0, nH):

                    imin=distances.index(min(distances))
                    distances.pop(imin)
                    hydrogens.append(indices.pop(imin))

        
                # Add molecule
                self.H2Oindx.append([iO] + hydrogens)


        return




    def EulerMat(self, alpha, beta, gamma):
        #
        # Euler rotation matrix defined as (Sakurai eq (3.3.18) )
        #
        # R(alpha,beta,gamma) = R_{z'}(gamma) * R_{y'}(beta) * R_{z}(alpha)
        #                     = R_{z}(alpha) * R_{y}(beta) * R_{z}(gamma) 
        #
        # INPUT        self.EulerAll(0.0, -theta, -phi)

        #    alpha, beta, gamma = Euler angles in radians
        #    
        # OUTPUT
        #    R
        
        R1 = np.array( [                 
                [np.cos(alpha), -np.sin(alpha), 0.0],
                [np.sin(alpha),  np.cos(alpha), 0.0],
                [0.0,            0.0,           1.0] 
                ] )
        
        R2 = np.array( [ 
                [np.cos(beta), 0.0, np.sin(beta)],
                [0.0,          1.0,           0.0],
                [-np.sin(beta), 0.0,  np.cos(beta)] 
                ] )

        R3 = np.array( [ 
                [np.cos(gamma), -np.sin(gamma), 0.0],
                [np.sin(gamma),  np.cos(gamma), 0.0],
                [0.0,            0.0,           1.0] 
                ] )
        

        # v' = R*v = R1*R2*R3*v is the rotated vector

        # Full rotation matrix
        R = np.dot(R1, np.dot(R2, R3) )
              
        return R




    def Euler(self, alpha, beta, gamma, vect):
        #
        # Euler rotation of vect defined as (Sakurai eq (3.3.18) )
        #
        # R(alpha,beta,gamma) = R_{z'}(gamma) * R_{y'}(beta) * R_{z}(alpha)
        #                     = R_{z}(alpha) * R_{y}(beta) * R_{z}(gamma) 
        #
        # INPUT        self.EulerAll(0.0, -theta, -phi)

        #    alpha, beta, gamma = Euler angles in radians
        #    vect
        #    
        # OUTPUT
        #    vect
        
        R1 = np.array( [                 
                [np.cos(alpha), -np.sin(alpha), 0.0],
                [np.sin(alpha),  np.cos(alpha), 0.0],
                [0.0,            0.0,           1.0] 
                ] )
        
        R2 = np.array( [ 
                [np.cos(beta), 0.0, np.sin(beta)],
                [0.0,          1.0,           0.0],
                [-np.sin(beta), 0.0,  np.cos(beta)] 
                ] )

        R3 = np.array( [ 
                [np.cos(gamma), -np.sin(gamma), 0.0],
                [np.sin(gamma),  np.cos(gamma), 0.0],
                [0.0,            0.0,           1.0] 
                ] )
        

        # v' = R*v = R1*R2*R3*v is the rotated vector

        # Full rotation matrix
        R = np.dot(R1, np.dot(R2, R3) )
      
        #aux = R1.dot( R2.dot( R3.dot(vect.T) ) )
        aux= np.dot(R,vect.T)
        vect = aux
        
        
        # http://stackoverflow.com/questions/986006/python-how-do-i-pass-a-variable-by-reference
        # Apparently Python passes arguments by assignment (discussion under link above). If we want 
        # to get out a modified quantity, the "return" command can be used to do so. Then the function
        # needs to be called by a statement of the type result=function(input)
        
        return vect




    def EulerDeg(self, alpha, beta, gamma, vect):
        #
        # Euler rotation of vect defined as (Sakurai eq (3.3.18) )
        #
        # R(alpha,beta,gamma) = R_{z'}(gamma) * R_{y'}(beta) * R_{z}(alpha)
        #                     = R_{z}(alpha) * R_{y}(beta) * R_{z}(gamma) 
        #
        # INPUT
        #    alpha, beta, gamma = Euler angles in degrees
        #    vect
        #    
        # OUTPUT
        #    vect
        
        R1 = np.array( [                 
                [np.cos(np.radians(alpha)), -np.sin(np.radians(alpha)), 0.0],
                [np.sin(np.radians(alpha)),  np.cos(np.radians(alpha)), 0.0],
                [0.0, 0.0, 1.0] ] )
        
        R2 = np.array( [ 
                [np.cos(np.radians(beta)), 0.0, np.sin(np.radians(beta))],
                [0.0,                      1.0, 0.0],
                [-np.sin(np.radians(beta)), 0.0, np.cos(np.radians(beta))] ] )

        R3 = np.array( [ 
                [np.cos(np.radians(gamma)), -np.sin(np.radians(gamma)),0.0],
                [np.sin(np.radians(gamma)),  np.cos(np.radians(gamma)), 0.0],
                [0.0, 0.0, 1.0] ] )
        


        # aux = R1*R2*R3*v is the rotated vector      
        aux = np.dot(R1, np.dot(R2, np.dot(R3,vect.T) ) )
        vect = aux
        
        # http://stackoverflow.com/questions/986006/python-how-do-i-pass-a-variable-by-reference
        # Apparently Python passes arguments by assignment (discussion under link above). If we want 
        # to get out a modified quantity, the "return" command can be used to do so. Then the function
        # needs to be called by a statement of the type result=function(input)
        
        return vect


    def EulerAll(self, alpha, beta, gamma):
        #
        # Euler rotation on all atoms
        #
        for iat in range(0, self.nat):
            aux=self.Euler(alpha, beta, gamma, self.x[iat])
            self.x[iat]=aux
        return



    def R(self, alpha, beta, gamma):
        #
        # Rotate entire system by R(alpha, beta, gamma)
        #
        # CAUTION: most probably the system wants to be translated to
        #          a reference origin before performing the rotation.
        #
        for iat in range(0, self.nat):
            self.x[iat] = Euler(alpha, beta, gamma, self.x[iat])

        return




    def ab(self, A, Ainv, a, b):
        #
        # Calculates the vector b-a with PBC
        #

        # Minimum image conversion                                                                                                              
        ab=b-a
        sab=np.dot(Ainv, ab)
        sab[:]=sab[:] - np.round(sab[:])
        ab=np.dot(A,sab)
        
        return ab




    def angle(self, v1, v2):
        #
        # Angle in radians between vectors v1, v2.
        # v1 and v2 must be numpy.array
        # 

        zero=np.array([0.0, 0.0, 0.0])

        v1norm=math.sqrt(np.dot(v1, v1))
        v2norm=math.sqrt(np.dot(v2, v2))
        v1dotv2=np.dot(v1, v2)


        alpha = math.acos(v1dotv2/(v1norm*v2norm))

        # Bring angle to [0, 2.0*pi[
        if (alpha < 0.0):
            alpha = alpha + 2.0*math.pi
           
        return alpha




    def angle_PBC(self, A, Ainv, x1, x0, x2):
        #
        # Angle in radians between points x1, x0, x2
        #   
        #  x1        x1 is one end
        #    \
        #     \
        #      \
        #       x0   x0 is the vertex
        #      /
        #     /
        #    /
        #  x2        x2 is the other end
        #
        #
        # 
        # x1, v and v2 must be numpy.array
        # 
        
        zero=np.array([0.0, 0.0, 0.0])
        
        v1=x1-x0
        v2=x2-x0
        
        
        
        # Minimum image conversion
        s1=np.dot(Ainv, v1)
        s1[:]=s1[:] - np.round(s1[:])
        v1=np.dot(A,s1)
        
        s2=np.dot(Ainv, v2)
        s2[:]=s2[:] - np.round(s2[:])
        v2=np.dot(A,s2)


        # Calculate angle
        v1norm=math.sqrt(np.dot(v1, v1))
        v2norm=math.sqrt(np.dot(v2, v2))
        v1dotv2=np.dot(v1, v2)

        ###print "DEBUG-- v1=", v1, "  v2=", v2
        alpha = math.acos(v1dotv2/(v1norm*v2norm))


        # Bring angle to [0, 2.0*pi[
        if (alpha < 0.0):
            alpha = alpha + 2.0*math.pi
           
        return alpha





    def printsnap(self, comment):
        #
        # Print snapshot in xyz format
        #
        print "\t", self.nat
        print comment
        for iat in range(0, self.nat):
            print self.tx[iat][0], "\t", self.x[iat][0], "\t", self.x[iat][1], "\t", self.x[iat][2]

        return




    def printsnapmod(self, u, comment):
        #
        # Print snapshot in xyz format plus an additional point u
        #
        print "\t", self.nat+1
        print comment
        for iat in range(0, self.nat):
            print self.tx[iat][0], "\t", self.x[iat][0], "\t", self.x[iat][1], "\t", self.x[iat][2]

        print "He", "\t", u[0], "\t", u[1], "\t", u[2]
            

        return

        


    def normalize(self, v):
        #
        # Normalize R^3 vector v
        #

        zero=np.array([0.0, 0.0, 0.0])

        # return normalized normal vector
        norm=self.d(zero, v)
        
        return v/norm 



    def bisector(self, v1, v2):
        #
        # Unit length bisector
        # v1 and v2 must be numpy.array
        #

        v12=v1+v2

        #zero=np.array([0.0, 0.0, 0.0])
        #return v12/self.d(zero, v12)


        
        return self.normalize(v12)
    



    def normal(self, v1, v2):
        #
        # Unit length normal vector.
        # v1 and v2 must be numpy.array
        #

        # Calculate cross product
        nv=np.cross(v1,v2)

        # Check that vectors were not collinear
        zero=np.array([0.0, 0.0, 0.0])
        if (self.d(zero,nv) < 1.0E-06):
            print " "
            print "WARNING!! Vectors are collinear and normal vector cannot be found. Do not trust results beyond this point."
            print " "

        # return normalized normal vector
        return self.normalize(nv)
        
    

    def phitheta(self, a):
        #
        # Given a unit vector a, find its spherical 
        # coordinates (phi, theta) in radians.
        # a must be numpy.array
        #

        
        # Checks on unit vector a
        
        tol=1.0E-12 # tolerance

        # make y component positive -- to avoid problems with arctan2()
        #if(a[1] < 0.0):
        #    a=-1.0*a

        if(abs(a[0]) < tol and abs(a[1]) < tol):
            phi=0.0
            theta=0.0
        
        else:
            phi = np.arctan2(a[1],a[0])
            theta = np.arccos(a[2])
        
        # Bring phi to [0, 2.0*pi[
            if (phi < 0.0):
                phi = phi + 2.0*math.pi
            if (theta >= 2.0*math.pi):
                theta = theta - 2.0*math.pi
        
        
        # Bring theta to [0, pi[
            if (theta < 0.0):
                theta = theta + 2.0*math.pi
        #   if (theta >= math.pi):
        #       theta = theta - math.pi
            
        
        return phi, theta




    def CaO(self, molindex):
        #
        # CaO stands for "Center and Orient"
        #
        # Center and orient system with respect
        # to the molecule labeled by molindex.
        #
        #
        # STEPS
        #
        # (1): Find indices for molecule at hand
        # (2): Bring oxygen to (0,0,0)
        # (3): Rotate so that the molecule bisector
        #         is along x and the normal along z
        #

        
       
        # (1)
        
        if not self.H2Oindx:
            self.FindH2Os
        
        iO=self.H2Oindx[molindex][0]
        iH1=self.H2Oindx[molindex][1]
        iH2=self.H2Oindx[molindex][2]
        
        
        # (2) 
        self.T(-1.0*self.x[iO]) # Translate by -x_O
        

        # (3): Rotate system by Ry[-theta]*Rz[-phi]
        u = self.normalize(self.x[iH1] - self.x[iO])
        v = self.normalize(self.x[iH2] - self.x[iO])
        

        normal=self.normal(u, v)
        phi, theta = self.phitheta(normal)
        
        #xaux=math.cos(phi)*math.sin(theta)
        #yaux=math.sin(phi)*math.sin(theta)
        #zaux=math.cos(theta)
        #aux= normal - np.array([xaux, yaux, zaux])

        self.EulerAll(0.0, -theta, -phi)



        bisector=self.bisector(self.x[iH1], self.x[iH2])
        
        alpha = self.angle(bisector, np.array([1.0, 0.0, 0.0] ) )

        
        # Consider all possible rotations
        if ( bisector[1] > 0.0):
            self.EulerAll(0.0, 0.0, -alpha)
        elif ( bisector[1] < 0.0):
            self.EulerAll(0.0, 0.0, alpha)


        return
        




    def xymirror(self):
        # Apply a mirror reflection about xy plane
        # to all atomic coordinates
       
        for iat in range(0, self.nat):
            self.x[iat][2] = -1.0*self.x[iat][2]

        return

            

    def molabovexy(self, imol):
        # If Oxygen of molecule with index imol is below the xy plane,
        # apply reflection about xy plane so that it lays above.
        
        iO=self.H2Oindx[imol][0]
        if(self.x[iO][2] < 0.0):
            self.xymirror()

        return


    
    def geteulerangles(self):
        # Get Euler angles of a selected molecule

        return


    


    
    def dimercoords(self, A, Ainv, imol1, imol2):
        # Get dimer internal coordinates
        #
        # rOH11, rOH12, rOH21, rOH22 - O-H distances (Ang)
        # HOH1, HOH2 - HOH angles (rad)
        # rOO, phi, theta - spherical coordinates of O2 
        # alpha, beta, gamma - Euler angles of molecule 2
        # 
        
        
        if (imol1 == imol2):
            
            print "ERROR: imol1 == imol2 in dimercoords(imol1, imol2). EXIT."
            exit()

        else:

            # intramolecular distances
            #rOH11=self.d(self.x[self.H2Oindx[imol1][0]], self.x[self.H2Oindx[imol1][1]])
            #rOH12=self.d(self.x[self.H2Oindx[imol1][0]], self.x[self.H2Oindx[imol1][2]])
            #rOH21=self.d(self.x[self.H2Oindx[imol2][0]], self.x[self.H2Oindx[imol2][1]])
            #rOH22=self.d(self.x[self.H2Oindx[imol2][0]], self.x[self.H2Oindx[imol2][2]])

            

            i11=3*imol1+1
            i12=3*imol1+2
            i21=3*imol2+1
            i22=3*imol2+2


            
            #rOH11=self.d_PBC(A, Ainv, self.x[self.H2Oindx[imol1][0]], self.x[i11] )
            #rOH12=self.d_PBC(A, Ainv, self.x[self.H2Oindx[imol1][0]], self.x[i12] )
            #rOH21=self.d_PBC(A, Ainv, self.x[self.H2Oindx[imol2][0]], self.x[i21] )
            #rOH22=self.d_PBC(A, Ainv, self.x[self.H2Oindx[imol2][0]], self.x[i22] )
            rOH11=d_PBC(A, Ainv, self.x[self.H2Oindx[imol1][0]], self.x[i11] )
            rOH12=d_PBC(A, Ainv, self.x[self.H2Oindx[imol1][0]], self.x[i12] )
            rOH21=d_PBC(A, Ainv, self.x[self.H2Oindx[imol2][0]], self.x[i21] )
            rOH22=d_PBC(A, Ainv, self.x[self.H2Oindx[imol2][0]], self.x[i22] )
            

            # intramolecular angles
            HOH1=self.angle_PBC(A, Ainv, self.x[i11] ,self.x[self.H2Oindx[imol1][0]] , self.x[i12]) 
            HOH2=self.angle_PBC(A, Ainv, self.x[i21] ,self.x[self.H2Oindx[imol2][0]] , self.x[i22]) 



            # O coordinates 
            # distance
            #rOO=self.d_PBC(A, Ainv, self.x[self.H2Oindx[imol1][0]], self.x[self.H2Oindx[imol2][0]])
            rOO=d_PBC(A, Ainv, self.x[self.H2Oindx[imol1][0]], self.x[self.H2Oindx[imol2][0]])

            # angles -- need to apply PBCs to correclty get direction vector
            u=self.x[self.H2Oindx[imol2][0]] - self.x[self.H2Oindx[imol1][0]]
            us=np.dot(Ainv, u)
            us[:]=us[:] - np.round(us[:])
            u=np.dot(A,us)
            u=self.normalize(u)
            
            phi, theta = self.phitheta(u)
            
            
            # Euler angles of secondary molecule
        
            # 1) Build rotation matrix
            
            iO=self.H2Oindx[imol2][0]
            iH1=self.H2Oindx[imol2][1]
            iH2=self.H2Oindx[imol2][2]


            # Atomic coordinates in system centered at secondary molecule        
            xH1aux = self.normalize(self.x[iH1] - self.x[iO])
            xH2aux = self.normalize(self.x[iH2] - self.x[iO])


            normal=self.normal(xH1aux,xH2aux)


            phiaux, thetaaux = self.phitheta(normal)
            R=self.EulerMat(0.0, -thetaaux, -phiaux)  # This matrix gives 2 out of 3 Euler angles      
        
            # Rotate normal of secondary molecule is along z
            xH1aux = np.dot(R,xH1aux)
            xH2aux = np.dot(R,xH2aux)
            
            # Calculate third Euler angle alpha
            bisector=self.bisector(xH1aux, xH2aux)
            alpha = self.angle(bisector, np.array([1.0, 0.0, 0.0] ) )
            
            # Euler angles are alpha, beta, gamma
            gamma=phiaux
            beta=thetaaux

            return [rOO, phi, theta, rOH11, rOH12, rOH21, rOH22, HOH1, HOH2, alpha, beta, gamma]
        


                
    
    def swap(self, i,j):
        #
        # Swap atoms i and j in the cluster
        #
        
        natpmol = self.nat/self.nmol
        
        molindxi = i // natpmol
        atindxi = i - molindxi * natpmol

        molindxj = j // natpmol
        atindxj = j - molindxj * natpmol

        self.H2Oindx[molindxi][atindxi], self.H2Oindx[molindxj][atindxj] = self.H2Oindx[molindxj][atindxj], self.H2Oindx[molindxi][atindxi] 
        self.tx[i], self.tx[j] = self.tx[j], self.tx[i]
        self.x[i], self.x[j] = self.x[j], self.x[i]
        
        return
        



    def findHB(self):
        #
        # Search for smallest INTERmolecular O-H distance
        # (corresponding to the H-bond if it exists) between 
        # molecules 1 and 2. 
        #
        # OUTPUT:
        # - centralmol labels the donor molecule (can take values 0,1) 
        # - secondarymol labels the acceptor molecule (can take values 0,1) 
        # - primaryh labels the donated H (can take values 1, 2, 4 or 5)
        # 


        # Find indices to all 6 atoms
        iO1=self.H2Oindx[0][0]
        iH11=self.H2Oindx[0][1]
        iH12=self.H2Oindx[0][2]
        iO2=self.H2Oindx[1][0]
        iH21=self.H2Oindx[1][1]
        iH22=self.H2Oindx[1][2]
        

        # Get all 4 INTERMOLECULAR O-H distances
        d1=self.d(self.x[iO1],self.x[iH21])
        d2=self.d(self.x[iO1],self.x[iH22])
        d3=self.d(self.x[iO2],self.x[iH11])
        d4=self.d(self.x[iO2],self.x[iH12])
        alld=[d1,d2,d3,d4]    

        
        # Minimum intramolecular O-H distance determines H bond.
        # Find it, assign corresponding O to secondary molecule 
        # (acceptor) and H to primary (donor) and then place 
        # set index values for molecules and Hydrogens


        minindx=alld.index(min(alld)) 

        if (minindx==0):
            centralmol=1
            secondarymol=0
            primaryh=iH21
            secondaryh=iH22

            
        elif(minindx==1):
            centralmol=1
            secondarymol=0
            primaryh=iH22
            secondaryh=iH21


        elif(minindx==2):
            centralmol=0
            secondarymol=1
            primaryh=iH11
            secondaryh=iH12

            
        elif(minindx==3):
            centralmol=0
            secondarymol=1
            primaryh=iH12
            secondaryh=iH11
                
 
            
        return centralmol, secondarymol, primaryh, secondaryh





    def isthisHB(self, imol1, imol2, HBdef, A, Ainv):
        # ASC -- HOMEBASE
        # Check if molecules imol1 and imol2 are forming an H-bond.
        # The bond is sought to be of the form
        #     O1-H1---O2
        # where imol1 is the donor and imol2 is the acceptor
        #
        #
        #
        #          HBdef values
        #
        # HBdef == 1 (FROM Corsetti et. al, JCP 139, 194502 (2013), Appendix B)
        # Needs to satisfy:
        #   (i) The intermolecular distance rOO < r_cut , where r_cut is 
        #       usually chosen as the position of OO OO the first minimum 
        #       in the O-O RDF (==3.5 Ang)
        #  (ii) The angle Oa-Od-Hd < theta^cut_OOH == 30deg
        #
        #
        # OUTPUT:
        #   HB ==  1 if imol1 and imol2 are H-bonded and imol1 is DONOR
        #   HB == -1 if imol1 and imol2 are H-bonded and imol1 is ACCEPTOR
        #   HB ==  0 if imol1 and imol2 are NOT H-bonded
        #

        # Check molecules are different
        if (imol1 == imol2):
            return 0




        # Initalize to False
        HB=0
        
        # Find indices to all 6 participating atoms
        iO1 =self.H2Oindx[imol1][0]
        iH11=self.H2Oindx[imol1][1]
        iH12=self.H2Oindx[imol1][2]
        iO2 =self.H2Oindx[imol2][0]
        iH21=self.H2Oindx[imol2][1]
        iH22=self.H2Oindx[imol2][2]

        
        if (HBdef==1):
            
            r_cut    = 3.5                  # (Ang)
            theta_cut= 30.0*self.deg2rad    # (rad)


            # Check (i) and (ii)
            if (self.D[iO1][iO2] <= r_cut): # Distance requirement passed.

                
                # Now establish possible H-bond and check angles.

                # Get all 4 INTERMOLECULAR O-H distances
                d121=self.D[iO1][iH21]
                d122=self.D[iO1][iH22]
                d211=self.D[iO2][iH11]
                d212=self.D[iO2][iH12]
                alld=[d121,d122,d211,d212]
                
                # Minimum intramolecular O-H distance determines possible 
                # H bond. Find it, assign corresponding O to acceptor molecule 
                # and H to donor molecule and then place set index values 
                # for molecules and Hydrogens.
                minindx=alld.index(min(alld)) 
                
                if (minindx==0):
                    donormol   = imol2
                    acceptormol= imol1
                    donorO     = iO2
                    acceptorO  = iO1
                    donorH     = iH21
                    passiveH   = iH22
            
                elif(minindx==1):
                    donormol   = imol2
                    acceptormol= imol1
                    donorO     = iO2
                    acceptorO  = iO1
                    donorH     = iH22
                    passiveH   = iH21
                    
                elif(minindx==2):
                    donormol   = imol1
                    acceptormol= imol2
                    donorO     = iO1
                    acceptorO  = iO2
                    donorH     = iH11
                    passiveH   = iH12
                    
                elif(minindx==3):
                    donormol    = imol1
                    acceptormol = imol2
                    donorO      = iO1
                    acceptorO   = iO2
                    donorH      = iH12
                    passiveH    = iH11
            

                # Now the possible H-bond is determined. Check angle Oa-Od-Hd angle.
                theta=self.angle_PBC(A, Ainv, self.x[acceptorO], self.x[donorO], self.x[donorH])
                
                if(theta <= theta_cut): # HB 
                    if  (imol1 == donormol):
                        HB = 1
                    elif(imol1 == acceptormol):
                        HB =-1
                else:
                    HB = 0
            
            else:
                HB = 0

        else:
            print "HBtype not recognized. Exiting..."
            exit()

            
        return HB



    def donHfirst(self, primaryh, secondaryh):
        # 
        # Make the donated Hydrogen be the first in the primary 
        # molecule 
        #
        # Call this function after calling findHB()
        #
        
        if (primaryh == 2 or primaryh == 5):

            self.swap(primaryh, secondaryh)
            primaryh, secondaryh = secondaryh, primaryh


        return primaryh, secondaryh

        


    def findHBandsort(self):
        #
        # Wrap functions findHB and donHfirst into a single one
        #
        centralmol, secondarymol, primaryh, secondaryh = self.findHB()
        primaryh, secondaryh = self.donHfirst(primaryh, secondaryh)

        return centralmol, secondarymol, primaryh, secondaryh
 




    def HBcoords(self, iD, iA, iH):
        # Given a donor O, an acceptor O and a 
        # donated H, find the relative coordinates
        #
        # nu = d(Od, Hd) - d(Oa, Hd)
        # mu = d(Od, Hd) + d(Oa, Hd)
        # OHO = angle (Od, H, Oa)
        #
        #
        # iD and iA are molecule indices, while 
        # iH is an atom index
        # 
        
        
        # Find Oxygen indices from molecule indices
        iDO = 3*iD ### + 0 # because O is the first atom in the triplet (O,H1,H2)
        iAO = 3*iA ### + 0

    
        # Get coordinates
        xOD = self.x[iDO]
        xOA = self.x[iAO]
        xH  = self.x[iH]
        
        dDH = self.d(xOD, xH)
        dAH = self.d(xOA, xH)
    
        
        v1 = np.array(xOD) - np.array(xH)
        v2 = np.array(xOA) - np.array(xH)

        nu = dDH - dAH
        mu = dDH + dAH
        #dODH = (mu+nu)/2.0
        #dOAH = (mu-nu)/2.0
        OHO = self.angle(v1, v2)
        

        return nu, mu, dDH, dAH, OHO




    def TOP(self, icentral, i1, i2, i3, i4):
        #
        # Tetrahedral Order Parameter (TOP)
        # q=1 - 3/8 \sum_{j=1}^3\sum_{i=1}^4 (cos \psi_{ij} + 1/3)^2
        #
        # icentral: Oxygen index w.r.t. which the TOP is to be calculated
        # i1,i2,i3,i4: indices of the 4 closest Oxygens to the central one
        #  
        #
        #
        
        onethird=1.0/3.0
        
        iext=[i1, i2, i3, i4]
        
        q=0.0
        
        # iO1 and iO2 run from 1 to 4
        # 
        iO1=0
        while iO1 < 3:
            
            iO2=iO1+1
            while iO2 < 4:
                
                q = q + (math.cos( self.angle(self.x[iext[iO1]] - self.x[icentral],x[iext[iO2]] - self.x[icentral])) + onethird)**2
                
                iO2+=1

            iO1+=1


        return 1.0 - 3.0/8.0 * q



    def getD(self, A, Ainv):
        #
        # Calculate D, the matrix of interatomic
        # distances. This is a symmetric matrix with
        # zero diagonal so we only need to calculate 
        # the distances for i2<i1
        #
        # Initialize matrix
        if (np.all(self.D)==False):
            self.D=np.zeros([self.nat,self.nat] )
        

        # Fill lower triangle of matrix
        for i1 in range(0,self.nat):
            i2=0
            while i2 < i1:
                if (i2==i1):
                    pass
                else:
                    self.D[i1][i2]=self.d_PBC(A, Ainv, i1, i2)
                
                i2+=1
        
        
        # Fill upper triangle of (symmetric, zero diagonal) matrix
        self.D=self.D + self.D.T
        
        return



    def OOdists(self, A, Ainv):
        #
        # Calculate matrix of O-O distances.
        #


        # If D has not been created yet
        if (np.size(self.D) != self.nat**2):

            if (np.size(self.rOO) != len(self.Oindx)**2 ):
                self.rOO=np.zeros( [len(self.Oindx), len(self.Oindx)] )
                

        
            # Fill lower triangle of matrix
            for iO1 in range(0,len(self.Oindx)):
                iO2=0
                while iO2 < iO1:
                    if (iO2==iO1):
                        pass
                    else:
                        self.rOO[iO1][iO2]=self.d_PBC(A, Ainv, self.Oindx[iO1], self.Oindx[iO2])
                
                    iO2+=1


            # Fill upper triangle of (symmetric, zero diagonal) matrix
            self.rOO=self.rOO + self.rOO.T
        

        # If D already exists
        else:
            # Take from D the O-O submatrix 
            self.rOO=self.D[self.Oindx,:][:,self.Oindx]
        
        return





    def OTO(self, A, Ainv, iO0, iO1to4):


        aux=[]

        # iO1 and iO2 run from 1 to 3 and 2 to 4 respectively

        iO1=0
        while (iO1 < 3):
            iO2=iO1+1
            while (iO2 < 4):
                aux.append(self.angle_PBC(A, Ainv,
                                            self.x[self.Oindx[iO1to4[iO1]]], self.x[iO0], self.x[self.Oindx[iO1to4[iO2]]]))

                iO2+=1

            iO1+=1

        aux=np.array(aux)
        aux=(1.0/3.0 + np.cos(aux))**2
        

        return 1.0 - (3.0/8.0)*aux.sum()



    def TTO(self, d1234):

        # iO1to4 are the distances of the closest 4 
        # oxygens to the central oxygen
        
        
        
        mean=np.sum(d1234)/4.0

        aux=d1234-mean

        return 1.0 - 1.0/(12.0*mean**2) * np.dot(aux,aux) 



    def getTOPs(self, A, Ainv, imolcentral):
        #
        # Wrapper to obtain Tetrahedral Order Parameters
        # - Orientational Tetrahedral Order (TOP)
        # - Translational Tetrahedral Order (TTP)
        #
        #
        # - A is the cell-->cartesian basis change matrix
        # - Ainv is the inverse of A
        # - imolcentral is the molecule index of the molecule
        #   in the cluster for which TOPs are being calculated
        #

        # Calculate O-O distances if needed
        if (np.size(self.rOO) != len(self.Oindx)**2 ):
            self.OOdists(A, Ainv)

    
        # calculate Orientational Tetrahedral Order parameter
        # self.rOO[imolcentral][:].argsort()[1:5] gets the
        # indices of the 4 smallest O-O distances
        q=self.OTO(A, Ainv, self.Oindx[imolcentral], self.rOO[imolcentral][:].argsort()[1:5] )



        # calculate Translational Tetrahedral Order parameter
        # np.sort(self.rOO[imolcentral][:])[1:5] takes the 4
        # smallest O-O distances for the given molecule
        Sk=self.TTO(np.sort(np.copy(self.rOO[imolcentral][:]))[1:5])

        return q, Sk





    def LSI(self, A, Ainv, imolcentral, rcutoff=3.7):
        # Calculate local structure index



        # Calculate O-O distances if needed                                                                                                              
        if (np.size(self.rOO) != len(self.Oindx)**2 ):
            self.OOdists(A, Ainv)



        # Select distances and sort them in place
        dists=np.copy(self.rOO[imolcentral][:])   # np.copy is important. Otherwise sort
        dists.sort()                              # command on this line will sort the
        dists=np.delete(dists,0)                  # original array!!



        # lesser than cutoff?
        ltc= dists < rcutoff

        # First molecule beyond cutoff radius needs to be included 
        nt=sum(ltc)
        ltc[nt]=True


        # Eliminate unnecessary array elements
        dists=dists[ltc]
        
        
        delta=dists[1:len(dists)] - dists[0:len(dists)-1]
        mean=np.sum(delta)/len(delta)

        return np.sum((delta-mean)**2)/len(delta)
        



    


#################################
# End of class Cluster
#################################



#################################
# Program functions
#################################


def file_len(fname):
    # From http://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
    # Count number of lines in file
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1        



def readANI(filename):
    #
    # Read ANI file. 
    # Returns a list with the number of atoms,
    # a list with all comment lines and a 
    # list of lists with all atomic positions.
    #
    # 
                                    
    nat=[]
    comment=[]
    txyz=[]


    numlines=file_len(filename)
    
    
    # Open file and read contents
    
    f=open(filename, 'r')
    
    txyzaux=[]
    j=1 # j takes values from 1 to 2+nat for each cluster
    ilin=1
    while ilin <= numlines:

        # j is an auxiliary index that helps keeping track
        # of what kind of read should be done.
        # j==1     : read number of atoms
        # j==2     : read comment line
        # 2<j<=nat : read atom type and coordinates
        row=f.readline()
        if (j==1):
            currentnat=int(row)
            nat.append(currentnat)
            j=j+1
                                               
        elif(j==2):
            currentcomment=str(row)
            comment.append(currentcomment)
            j=j+1
            
        else:
            txyzaux.append(row.split())
            if (j==currentnat+2):
                txyz.append(txyzaux)
                txyzaux=[]
                j=1
            else:
                j=j+1
            
        ilin = ilin + 1


        nsnap=len(nat)

    return nat, comment, txyz



def d(x1, x2):
    #
    # Cartesian distance between points x1 and x2.
    #
    return math.sqrt( (x1[0]-x2[0])**2 + (x1[1]-x2[1])**2 + (x1[2]-x2[2])**2 )



def d_PBC(A, Ainv, x1, x2):
    #
    # Cartesian distance between 2 points
    # in a system with periodic boundary
    # conditions (PBC)
    #
    # x1, x2: numpy arrays containing cartesian coordinates
    # A: numpy matrix containing cell vectors
    # Ainv: inverse of A
    #
    # Steps:
    #
    # 1) calculate vector linking x1 and x2, r12
    # 2) change to cell basis
    # 3) translate to mimimum image
    # 4) change to cartesian basis
    # 5) calculate euclidean distance
    # 
    
    r12=x1-x2
    s12=np.dot(Ainv, r12)
    s12[:]=s12[:] - np.round(s12[:])
    r12=np.dot(A,s12)
    return np.sqrt( np.dot(r12,r12) )




def printdimerx(x):
    # Print dimer coordinates
    # 
    # dimerx is a list containing
    # [rOO, phi, theta, rOH11, rOH12, rOH21, rOH22, HOH1, HOH2, alpha, beta, gamma]
    # or other coordinates
    #
    
    numcoords=len(x)
    
    entrylen=8 # lenght of each entry in output

    #print ["%0.6f" % i for i in x]    
    #[out.append(str(x[i])) for i in range(0, numcords)]
    
    out=""
    sep="\t"
    for i in range(0, numcoords):
        if (i==numcoords):
            sep=""
        out = out + str(x[i])[0:entrylen] + sep

    print out
    

    return 
