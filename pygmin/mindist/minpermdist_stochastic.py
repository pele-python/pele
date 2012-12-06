import numpy as np

from pygmin.mindist import CoMToOrigin, aa2xyz, alignRotation, findBestPermutation, getDistxyz
from pygmin.utils.rotations import random_aa, aa2mx, q2mx
from pygmin.mindist.mindistutils import getAlignRotation
from pygmin.mindist import ExactMatchCluster
from pygmin.mindist.distpot import MinPermDistPotential
import pygmin.defaults as defaults

__all__ = ["minPermDistStochastic"]

def applyRotation(mx, X1d):
    X = X1d.reshape([-1,3])
    X = np.dot(mx, X.transpose()).transpose()
    return X.reshape(-1)

def _optimizePermRot(X1, X2, niter, permlist, verbose=False, use_quench=True):
    if use_quench:
        pot = MinPermDistPotential(X1, X2.copy(), permlist=permlist)

    distbest = getDistxyz(X1, X2)
    mxbest = np.identity(3)
    X20 = X2.copy()
    for i in range(niter):
        #get and apply a random rotation
        aa = random_aa()
        if not use_quench:
            mx = aa2mx(aa)
            mxtot = mx
            #print "X2.shape", X2.shape
        else:
            #optimize the rotation using a permutationally invariand distance metric
            ret = defaults.quenchRoutine(aa, pot.getEnergyGradient, tol=0.01)
            aa1 = ret[0]
            mx1 = aa2mx(aa1)
            mxtot = mx1
        X2 = applyRotation(mxtot, X20)
        
        #optimize the permutations
        dist, X1, X2 = findBestPermutation(X1, X2, permlist)
        if verbose:
            print "dist", dist, "distbest", distbest
        #print "X2.shape", X2.shape
        
        #optimize the rotation
        dist, Q2 = getAlignRotation(X1, X2)
#        print "dist", dist, "Q2", Q2
        mx2 = q2mx(Q2)
        mxtot = np.dot(mx2, mxtot)
        
        if dist < distbest:
            distbest = dist
            mxbest = mxtot
    return distbest, mxbest
    
    

def minPermDistStochastic(X1, X2, niter=100, permlist=None, verbose=False, accuracy=0.01,
                      check_inversion=True, use_quench=False):
    """
    Minimize the distance between two clusters.  
    
    Parameters
    ----------
    X1, X2 : 
        the structures to align.  X2 will be aligned with X1, both
        the center of masses will be shifted to the origin
    niter : int
        the number of basinhopping iterations to perform
    permlist : a list of lists of atoms 
        A list of lists of atoms which are interchangable.
        e.g. if all the atoms are interchangable
        
            permlist = [range(natoms)]
        
        For a 50/50 binary mixture, 
        
            permlist = [range(1,natoms/2), range(natoms/2,natoms)]
    verbose : 
        whether to print status information
    accuracy : 
        accuracy for determining if the structures are identical
    check_inversion :
        if true, account for point inversion symmetry
    use_quench : 
        for each step of the iteration, minimize a permutationally invariant
        distance metric.  This slows the algorithm, but can potentially make
        it more accurate.

    Notes
    -----

    The following symmetries will be accounted for
    
        Translational symmetry

        Global rotational symmetry

        Permutational symmetry
        
        Point inversion symmetry

    
    This method should have the same outcome as minPermDistStochastic, but 
    uses a different method.  The algorithm here to find the best distance is
    
    for i in range(niter):    
        random_rotation(coords)
        findBestPermutation(coords)
        alignRotation(coords)
    """
    natoms = len(X1) / 3
    if permlist is None:
        permlist = [range(natoms)]

    X1init = X1
    X2init = X2
    X1 = np.copy(X1)
    X2 = np.copy(X2)

    #first check for exact match
    exactmatch = ExactMatchCluster(accuracy=accuracy, permlist=permlist)
    if exactmatch(X1, X2):
        #this is kind of cheating, I would prefer to return
        #X2 in best alignment and the actual (small) distance
        return 0.0, X1, X1.copy() 
    
    #bring center of mass of x1 and x2 to the origin
    #save the center of mass of X1 for later
    X1com = X1.reshape([-1,3]).sum(0) / natoms
    X1 = CoMToOrigin(X1)
    X2 = CoMToOrigin(X2)
    #print "X2.shape", X2.shape
    
    #find the best rotation stochastically
    X20 = X2.copy()
    distbest, mxbest = _optimizePermRot(X1, X2, niter, permlist, verbose=verbose, use_quench=use_quench)
    use_inversion = False
    if check_inversion:
        X20i = -X20.copy()
        X2 = X20i.copy()
        distbest1, mxbest1 = _optimizePermRot(X1, X2, niter, permlist, verbose=verbose, use_quench=use_quench)
        if distbest1 < distbest:
            if verbose:
                print "using inversion in minpermdist"
            use_inversion = True
            distbest = distbest1
            mxbest = mxbest1

    #now we know the best rotation
    if use_inversion: X20 = X20i
    X2 = applyRotation(mxbest, X20)
    dist, X1, X2 = findBestPermutation(X1, X2, permlist)
    dist, X2 = alignRotation(X1, X2)
    if dist > distbest+0.001:
        print "ERROR: minPermDistRanRot: dist is different from distbest %f %f" % (dist, distbest)
    if verbose:
        print "finaldist", dist, "distmin", distbest
    
    #add back in the center of mass of X1
    X1 = X1.reshape([-1,3])
    X2 = X2.reshape([-1,3])
    X1 += X1com
    X2 += X1com
    X1 = X1.reshape(-1)
    X2 = X2.reshape(-1)
    
    return dist, X1, X2


#
#
# below here only testing stuff
#
#

import unittest
from testmindist import TestMinDist
class TestMinPermDistStochastic_BLJ(TestMinDist):
    def setUp(self):
        from pygmin.potentials.ljpshiftfast import LJpshift as BLJ
        from pygmin import defaults
        
        self.natoms = 25
        self.ntypeA = int(self.natoms * .8)
        self.pot = BLJ(self.natoms, self.ntypeA)
        self.permlist = [range(self.ntypeA), range(self.ntypeA, self.natoms)]
        
        self.X1 = np.random.uniform(-1,1,[self.natoms*3])*(float(self.natoms))**(1./3)/2
        
        #run a quench so the structure is not crazy
        ret = defaults.quenchRoutine(self.X1, self.pot.getEnergyGradient, **defaults.quenchParams)
        self.X1 = ret[0]
        

    def testBLJ(self):
        import pygmin.defaults
        X1 = np.copy(self.X1)
        X2 = np.random.uniform(-1,1,[self.natoms*3])*(float(self.natoms))**(1./3)/2
        
        #run a quench so the structure is not crazy
        ret = pygmin.defaults.quenchRoutine(X2, self.pot.getEnergyGradient)
        X2 = ret[0]

        self.runtest(X1, X2, minPermDistStochastic)


    def testBLJ_isomer(self):
        """
        test with BLJ potential.  We have two classes of permutable atoms  
        
        test case where X2 is an isomer of X1.
        """
        import pygmin.utils.rotations as rot
        X1i = np.copy(self.X1)
        X1 = np.copy(self.X1)        
        X2 = np.copy(X1)
        
        #rotate X2 randomly
        aa = rot.random_aa()
        rot_mx = rot.aa2mx( aa )
        for j in range(self.natoms):
            i = 3*j
            X2[i:i+3] = np.dot( rot_mx, X1[i:i+3] )
        
        #permute X2
        import random, copy
        from pygmin.mindist.permutational_alignment import permuteArray
        for atomlist in self.permlist:
            perm = copy.copy(atomlist)
            random.shuffle( perm )
            X2 = permuteArray( X2, perm)

        X2i = np.copy(X2)
        
        #distreturned, X1, X2 = self.runtest(X1, X2)
        distreturned, X1, X2 = self.runtest(X1, X2, minPermDistStochastic)

        
        #it's an isomer, so the distance should be zero
        self.assertTrue( abs(distreturned) < 1e-14, "didn't find isomer: dist = %g" % (distreturned) )


def test(X1, X2, lj, atomtypes=["LA"], permlist = None, fname = "lj.xyz",
         minPermDist=minPermDistStochastic):
    import copy
    natoms = len(X1) / 3
    if permlist == None:
        permlist = [range(natoms)]
    
    X1i = copy.copy(X1)
    X2i = copy.copy(X2)
    
    printlist = []
    printlist.append((X2.copy(), "X2 initial"))
    printlist.append((X1.copy(), "X1 initial"))


    distinit = np.linalg.norm(X1-X2)
    print "distinit", distinit

    (dist, X1, X2) = minPermDist(X1,X2, permlist=permlist)
    distfinal = np.linalg.norm(X1-X2)
    print "dist returned    ", dist
    print "dist from coords ", distfinal
    print "initial energies (post quench)", lj.getEnergy(X1i), lj.getEnergy(X2i)
    print "final energies                ", lj.getEnergy(X1), lj.getEnergy(X2)

    printlist.append((X1.copy(), "X1 final"))
    printlist.append((X2.copy(), "X2 final"))


    import pygmin.printing.print_atoms_xyz as printxyz
    with open(fname, "w") as fout:
        for xyz, line2 in printlist:
            printxyz.printAtomsXYZ(fout, xyz, line2=line2 +" "+ str(lj.getEnergy(xyz)))

    

def test_binary_LJ(natoms = 12, **kwargs):
    import pygmin.defaults
    import pygmin.utils.rotations as rot
    quench = pygmin.defaults.quenchRoutine

    printlist = []
    
    ntypea = int(natoms*.8)
    from pygmin.potentials.ljpshift import LJpshift
    lj = LJpshift(natoms, ntypea)
    permlist = [range(ntypea), range(ntypea, natoms)]

    X1 = np.random.uniform(-1,1,[natoms*3])*(float(natoms))**(1./3)/2
    printlist.append( (X1.copy(), "very first"))
    #quench X1
    ret = quench( X1, lj.getEnergyGradient)
    X1 = ret[0]
    printlist.append((X1.copy(), "after quench"))

    X2 = np.random.uniform(-1,1,[natoms*3])*(float(natoms))**(1./3)
    #make X2 a rotation of X1
    print "testing with", natoms, "atoms,", ntypea, "type A atoms, with X2 a rotated and permuted isomer of X1"
    aa = rot.random_aa()
    rot_mx = rot.aa2mx( aa )
    for j in range(natoms):
        i = 3*j
        X2[i:i+3] = np.dot( rot_mx, X1[i:i+3] )
    printlist.append((X2.copy(), "x2 after rotation"))
    

    

    import random, copy
    from pygmin.mindist.permutational_alignment import permuteArray

    for atomlist in permlist:
        perm = copy.copy(atomlist)
        random.shuffle( perm )
        print perm
        X2 = permuteArray( X2, perm)
    printlist.append((X2.copy(), "x2 after permutation"))


    #X1 = np.array( [ 0., 0., 0., 1., 0., 0., 0., 0., 1.,] )
    #X2 = np.array( [ 0., 0., 0., 1., 0., 0., 0., 1., 0.,] )
    X1i = copy.copy(X1)
    X2i = copy.copy(X2)
    
    atomtypes = ["N" for i in range(ntypea)]
    for i in range(natoms-ntypea):
        atomtypes.append("O")
    
    print "******************************"
    print "testing binary LJ  ISOMER"
    print "******************************"
    test(X1, X2, lj, atomtypes=atomtypes, permlist = permlist, **kwargs)
    
    print "******************************"
    print "testing binary LJ  non isomer"
    print "******************************"
    X2 = np.random.uniform(-1,1,[natoms*3])*(float(natoms))**(1./3)
    ret = quench( X2, lj.getEnergyGradient)
    X2 = ret[0]
    test(X1, X2, lj, atomtypes=atomtypes, permlist=permlist, **kwargs)

    
        
        
def test_LJ(natoms = 12, **kwargs):
    from pygmin.potentials.lj import LJ
    import pygmin.defaults
    import pygmin.utils.rotations as rot
    from pygmin.mindist.permutational_alignment import permuteArray
    import random

    quench = pygmin.defaults.quenchRoutine
    lj = LJ()
    X1 = np.random.uniform(-1,1,[natoms*3])*(float(natoms))**(1./3)
    #quench X1
    ret = quench( X1, lj.getEnergyGradient)
    X1 = ret[0]
    X2 = np.random.uniform(-1,1,[natoms*3])*(float(natoms))**(1./3)
    #make X2 a rotation of X1
    print "testing with", natoms, "atoms, with X2 a rotated and permuted isomer of X1"
    aa = rot.random_aa()
    rot_mx = rot.aa2mx( aa )
    for j in range(natoms):
        i = 3*j
        X2[i:i+3] = np.dot( rot_mx, X1[i:i+3] )
    perm = range(natoms)
    random.shuffle( perm )
    print perm
    X2 = permuteArray( X2, perm)

    #X1 = np.array( [ 0., 0., 0., 1., 0., 0., 0., 0., 1.,] )
    #X2 = np.array( [ 0., 0., 0., 1., 0., 0., 0., 1., 0.,] )
    import copy
    X1i = copy.copy(X1)
    X2i = copy.copy(X2)
    
    print "******************************"
    print "testing normal LJ  ISOMER"
    print "******************************"
    test(X1, X2, lj, **kwargs)
    
    print "******************************"
    print "testing normal LJ  non isomer"
    print "******************************"
    X2 = np.random.uniform(-1,1,[natoms*3])*(float(natoms))**(1./3)
    ret = quench( X2, lj.getEnergyGradient)
    X2 = ret[0]
    test(X1, X2, lj, **kwargs)
    

    distinit = np.linalg.norm(X1-X2)
    print "distinit", distinit



if __name__ == "__main__":
    print "******************************"
    print "testing normal LJ"
    print "******************************"
    test_LJ(12)
    print ""
    print ""
    print "************************************"
    print "testing binary LJ with permute lists"
    print "************************************"
    test_binary_LJ(12)
    
    unittest.main()
    
    
