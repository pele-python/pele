import numpy as np
import copy
import rotations as rot
import potentials.potential as potential


def overlap_slow( XA, XB, L2, atomlist, nlist = [] ):
    E = 0.
    for i1 in atomlist:
        i = i1*3
        for j1 in atomlist:
            j = j1*3
            r2 = np.sum( (XB[i:i+3] - XA[j:j+3])**2 )
            E += np.exp(-r2/L2)
    return E

def overlap_gradient_slow( XA, XB, L2, atomlist, nlist = [] ):
    E = 0.
    grad = np.zeros(len(XA))
    iL2 = 1./L2
    for i1 in atomlist:
        i = i1*3
        for j1 in atomlist:
            j = j1*3
            dr = (XB[i:i+3] - XA[j:j+3])
            r2 = np.sum( dr**2 )
            de = np.exp(-r2/L2)
            E += de
            grad[ i:i+3 ] -= dr[:] * (de * 2.0 * iL2)
            #grad[ :j+3 ] += dr[:] * (de * 2.0 * iL2)
    return E, grad

def overlap_fast( XA, XB, L2, atomlist, nlist = [] ):
    xa = XA.reshape(XA.size/3,3)[atomlist,:]
    xb = XB.reshape(XB.size/3,3)[atomlist,:]
    #note: drmat = xa[:,np.newaxis] - xb[:] is an array of shape (natoms,natoms,3)
    #      drmat[i,j,:] == xa[i,:] - xb[j,:]
    return np.sum(np.exp(-np.sum((xa[:,np.newaxis] - xb[:])**2, axis=2)/L2))
   
class MinPermDistPotential(potential.potential):
    """
    Find the rotation (in angle-axis representation) which maximizes the
    permutation independent overlap between two structures.

    This potential defines the energy as the overlap, here defined as

    E = - sum_i sum_j exp( -|x_Ai - X_Bj|**2 / L**2 )

    There are other options for the overlap which are not implemented.  For
    instance something with a slower decay than exp(-R**2).
    """
    def __init__(self, XA, XB, L=0.2, permlist = []):
        self.XA0 = copy.copy(XA)
        self.XB0 = copy.copy(XB)
        self.XB = copy.copy(XB)
        self.nsites = len(XA)/3
        self.L2 = L*L
        self.permlist = permlist
        if len(self.permlist) == 0:
            self.permlist = [range( self.nsites)]
        try:
            from overlap import overlap, overlap_gradient
            self.overlap = overlap
            self.overlap_gradient = overlap_gradient
        except:
            #self.overlap = overlap_fast
            self.overlap = overlap_fast
            self.overlap_gradient = overlap_gradient_slow
            print "Using python energy calculation. Compile overlap.f90 to speed things up (a little)"

        self.setUpRigidBodies()

    def setUpRigidBodies(self):
        #set up one rigid body
        from potentials.rigid_bodies.molecule import Molecule
        self.rb = Molecule()
        for atomlist in self.permlist:
            #make a rigid body from atomlist
            for i in atomlist:
                self.rb.insert_site( 0, self.XB0[i*3:i*3+3] )

    def getEnergy_no_rigid_body(self, AA ):
        # Rotate XB0 according to angle axis AA
        rot_mx = rot.aa2mx( AA )
        for j in range(self.nsites):
            i = 3*j
            self.XB[i:i+3] = np.dot( rot_mx, self.XB0[i:i+3] )
        #return distance between XB and XA0
        E = 0.
        for atomlist in self.permlist:
            E -= self.overlap( self.XA0, self.XB, self.L2, atomlist, [self.nsites, len(atomlist)])
            #if E < -10.5:
            #   print E, atomlist, [self.nsites, len(atomlist)], len(self.XA0), len(self.XB)
        return E

    def getEnergy(self, AA):
        # Rotate rigid body according to angle axis AA
        self.XB = self.rb.getxyz( aa = AA )
        #get energy and gradient for each site
        E = 0.
        for atomlist in self.permlist:
            E -= self.overlap( self.XA0, self.XB, self.L2, atomlist, [self.nsites, len(atomlist)])
        return E


    def getEnergyGradient(self, AA):
        # Rotate XB0 according to angle axis AA
        # Rotate rigid body according to angle axis AA
        self.XB = self.rb.getxyz( aa = AA )
        #get energy and gradient for each site
        E = 0.
        grad = np.zeros(len(self.XB))
        for atomlist in self.permlist:
            de, dg = self.overlap_gradient( self.XA0, self.XB, self.L2, atomlist, [self.nsites, len(atomlist)])
            E -= de
            grad -= dg
        #convert site-gradients into angle -axis gradients
        #first calculate the rotation matrix and derivatives
        comgrad, aagrad = self.rb.getGradients(AA, grad, True)
        return E, aagrad


    def globalEnergyMin(self):
        """
        return the lowest energy theoretically possible.  This will happen if XB == XA
        """
        E = 0.
        for atomlist in self.permlist:
            E -= self.overlap( self.XA0, self.XA0, self.L2, atomlist, [self.nsites, len(atomlist)])
        return E


class RandomRotationTakeStep(object):
    def takeStep(self, aa):
        return random_rotation(aa)
    def updateStep(self, accept):
        pass

def random_rotation( aa):
    aanew = rot.random_aa()
    #print "aanew ", aanew
    aa[:] = aanew[:]
    #print "aa    ", aa
    

def test_distpot():
    #define two structures
    natoms = 12
    X1 = np.random.uniform(-1,1,[natoms*3])*(float(natoms))**(1./3)
    
    X2 = np.random.uniform(-1,1,[natoms*3])*(float(natoms))**(1./3)
    #make X2 a rotation of X1
    print "testing with", natoms, "atoms, with X2 a rotated and permuted isomer of X1"
    aa = rot.random_aa()
    rot_mx = rot.aa2mx( aa )
    for j in range(natoms):
        i = 3*j
        X2[i:i+3] = np.dot( rot_mx, X1[i:i+3] )
    #import random, mindistutils
    #perm = range(natoms)
    #random.shuffle( perm )
    #print perm
    #X2 = mindistutils.permuteArray( X2, perm)

    pot = MinPermDistPotential(X1, X2, L = .2)
    
    aa = rot.random_aa()
    e = pot.getEnergy(aa)
    print "energy", e
    de, dg = pot.getEnergyGradient(aa)
    print "energy from gradient", de, "diff", e-de
    
    den, dgn = pot.getEnergyGradientNumerical(aa)
    maxgrad= np.max( np.abs( dg ) )
    maxdiff = np.max( np.abs( dg - dgn ) )
    print "maximum gradient", maxgrad
    print "max difference in analytical vs numerical gradient", maxdiff


if __name__ == "__main__":
    test_distpot()
