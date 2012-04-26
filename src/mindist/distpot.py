import numpy as np
import copy
import rotations as rot
import potentials.potential as potential

def overlap_slow( XA, XB, L2, nsites ):
    E = 0.
    for i1 in range(nsites):
        i = i1*3
        for j1 in range(nsites):
            j = j1*3
            r2 = np.sum( (XB[i:i+3] - XA[j:j+3])**2 )
            E += np.exp(-r2/L2)
    return E

class MinPermDistPotential(potential.potential):
    """
    Find the rotation (in angle-axis representation) which maximizes the
    permutation independent overlap between two structures.

    This potential defines the energy as the overlap, here defined as

    E = - sum_i sum_j exp( -|x_Ai - X_Bj|**2 / L**2 )

    There are other options for the overlap which are not implimented.  For
    instance something with a slower decay than exp(-R**2).

    TODO: impliment analytical Gradients
    """
    def __init__(self, XA, XB, L=0.2):
        self.XA0 = copy.copy(XA)
        self.XB0 = copy.copy(XB)
        self.XB = copy.copy(XB)
        self.nsites = len(XA)/3
        self.L2 = L*L
        try:
            from overlap import overlap
            self.overlap = overlap
        except:
            self.overlap = overlap_slow
            print "Using slow energy calculation. Compile overlap.f90 to speed things up"
            

    def getEnergy(self, AA ):
        # Rotate XB0 according to angle axis AA
        rot_mx = rot.aa2mx( AA )
        for j in range(self.nsites):
            i = 3*j
            self.XB[i:i+3] = np.dot( rot_mx, self.XB0[i:i+3] )
        #return distance between XB and XA0
        E = -self.overlap( self.XA0, self.XB, self.L2, self.nsites)
        return E

    def globalEnergyMin(self):
        """
        return the lowest energy theoretically possible.  This will happen if XB == XA
        """
        E = -self.overlap( self.XA0, self.XA0, self.L2, self.nsites)
        return E



def random_rotation( aa):
    aanew = rot.random_aa()
    #print "aanew ", aanew
    aa[:] = aanew[:]
    #print "aa    ", aa
