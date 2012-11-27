import numpy as np

from mindistutils import CoMToOrigin, aa2xyz, alignRotation, findBestPermutation, getDistxyz
from pygmin.utils.rotations import random_aa, aa2mx, q2mx
from pygmin.mindist.mindistutils import getAlignRotation
from pygmin.mindist import ExactMatchCluster
from pygmin.mindist.distpot import MinPermDistPotential
import pygmin.defaults as defaults

__all__ = ["minPermDistRanRot"]

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
    
    

def minPermDistRanRot(X1, X2, niter=100, permlist=None, verbose=False, accuracy=0.01,
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

def testranrot_lj(natoms=38, **kwargs):
    from pygmin.potentials.lj import LJ
    from pygmin.optimize.quench import mylbfgs as quench
    from pygmin.mindist.minpermdist_stochastic import test, minPermDistStochastic
    import time

    lj = LJ()
    X1 = np.random.uniform(-1,1,[natoms*3])*(float(natoms))**(1./3)
    #quench X1
    ret = quench( X1, lj.getEnergyGradient)
    X1 = ret[0]
    X2 = np.random.uniform(-1,1,[natoms*3])*(float(natoms))**(1./3)
    #make X2 a rotation of X1
    print "testing with", natoms, "atoms, with X2 a rotated and permuted isomer of X1"
    aa = random_aa()
    rot_mx = aa2mx( aa )
    for j in range(natoms):
        i = 3*j
        X2[i:i+3] = np.dot( rot_mx, X1[i:i+3] )
    import random, mindistutils
    perm = range(natoms)
    random.shuffle( perm )
    print perm
    X2 = mindistutils.permuteArray( X2, perm)

    #X1 = np.array( [ 0., 0., 0., 1., 0., 0., 0., 0., 1.,] )
    #X2 = np.array( [ 0., 0., 0., 1., 0., 0., 0., 1., 0.,] )
    import copy
    X1i = copy.copy(X1)
    X2i = copy.copy(X2)
    
    print "******************************"
    print "testing normal LJ  ISOMER"
    print "******************************"
    print ""
    print "results from minPermDistRanRot"
    print ""
    X1, X2 = X1i.copy(), X2i.copy()    
    test(X1, X2, lj, minPermDist=minPermDistRanRot)
    print ""
    print "results from minPermDistStochastic"
    print ""
    X1, X2 = X1i.copy(), X2i.copy()
    test(X1, X2, lj, minPermDist=minPermDistStochastic)
    
    print "******************************"
    print "testing normal LJ  non isomer"
    print "******************************"
    X2 = np.random.uniform(-1,1,[natoms*3])*(float(natoms))**(1./3)
    ret = quench( X2, lj.getEnergyGradient)
    X2 = ret[0]
    X2i = X2.copy()
    
    print ""
    print "results from minPermDistRanRot"
    print ""
    X1, X2 = X1i.copy(), X2i.copy()  
    t1 = time.clock()  
    test(X1, X2, lj, minPermDist=minPermDistRanRot)
    t2 = time.clock()
    print "time elapsed", t2-t1
    print ""
    print "results from minPermDistStochastic"
    print ""
    X1, X2 = X1i.copy(), X2i.copy()
    t1 = time.clock()
    test(X1, X2, lj, minPermDist=minPermDistStochastic)
    t2 = time.clock()
    print "time elapsed", t2-t1
    

    #test(X1, X2, lj, minPermDist=minPermDistRanRot) 
    
if __name__ == "__main__":
    from pygmin.mindist.minpermdist_stochastic import test_LJ, test_binary_LJ
    print "******************************"
    print "testing normal LJ"
    print "******************************"
    test_LJ(12, minPermDist=minPermDistRanRot)
    print ""
    print ""
    print "************************************"
    print "testing binary LJ with permute lists"
    print "************************************"
    test_binary_LJ(12)
    
    print ""
    print ""
    print "testing ranrot"
    print ""
    print ""
    testranrot_lj()
#    unittest.main()
    
    
