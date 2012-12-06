import numpy as np
import itertools


__all__ = ["findBestPermutation"] 


_hungarian = None

if _hungarian is None:
    try:
        #use the hungarian package which is compiled
        #note, as of version 0.2.1 hungarian is buggy.
        # you can get a fixed version from 
        # https://github.com/js850/hungarian
        # hopefully the fixes will be incorporated into the official package soon
        import hungarian
        def _hungarian(cost):
            newind1 = hungarian.lap(cost)
            return [(i, j) for i,j in enumerate(newind1[0])]
        #print "hungari newind", newind
    except ImportError:
        pass

# if we didn't find hungarian algorithm before, try munkres
if _hungarian is None:
    try:
        #use the munkres package
        #convert cost matrix to a form used by munkres
        from munkres import Munkres
        
        def _hungarian(cost):
            matrix = cost.tolist()
            m = Munkres()
            return m.compute(matrix)
        #print "munkres newind", newind
    except ImportError:
        pass
    
if _hungarian is None:
    raise BaseException("No Hungarian algorithm implementation found!"
                    "Please install the hungarian or the munkres package")



def permuteArray(Xold, perm):
    #don't modify Xold
    Xnew = np.copy(Xold)
    permsorted = sorted(perm)
    for (iold, inew) in itertools.izip(permsorted, perm):
        #print iold, "->", inew
        Xnew[inew*3:inew*3+3] = Xold[iold*3:iold*3+3]

    return Xnew


def makeCostMatrix(X1, X2, atomlist):
    """
    return the cost matrix for use in the hungarian algorithm.
    
    the cost matrix is the distance matrix (squared) for all atoms in atomlist
    """
    atomlistnp = np.array(atomlist)
    X13 = X1[atomlistnp,:]
    X23 = X2[atomlistnp,:]
    cost = (((X13[np.newaxis,:] - X23[:,np.newaxis,:])**2).sum(2))
    return cost


def findBestPermutationList( X1, X2, atomlist = None, cost_function = None ):
    """
    For a given set of positions X1 and X2, find the best permutation of the
    atoms in X2.

    Use an implementation of the Hungarian Algorithm in the Python package
    index (PyPi) called munkres (another name for the algorithm).  The
    hungarian algorithm time scales as O(n^3), much faster than the O(n!) from
    looping through all permutations.

    http://en.wikipedia.org/wiki/Hungarian_algorithm
    http://pypi.python.org/pypi/munkres/1.0.5.2
    
    another package, hungarian, implements the same routine in comiled C
    http://pypi.python.org/pypi/hungarian/
    When I first downloaded this package I got segfaults.  The problem for me
    was casing an integer pointer as (npy_intp *).  I may add the corrected 
    version to pygmin at some point
    """
    nsites = len(X1) / 3

    if atomlist == None:
        atomlist = range(nsites)
    nperm = len(atomlist)

    #print "atomlist", atomlist

    #########################################
    # create the cost matrix
    # cost[j,i] = (X1(i,:) - X2(j,:))**2
    #########################################
#    cost = np.zeros( [nperm,nperm], np.float64)
#    for i in range(nperm):
#        atomi = atomlist[i]
#        for j in range(nperm):
#            atomj = atomlist[j]
#            R2 = np.sum( (X1[atomi*3:atomi*3+3] - X2[atomj*3:atomj*3+3])**2 )
#            cost[j,i] = R2

    X1 = X1.reshape([-1,3])
    X2 = X2.reshape([-1,3])
    cost = makeCostMatrix(X1, X2, atomlist)
    #cost = np.sqrt(cost)

    #########################################
    # run the hungarian algorithm
    #########################################
    newind = _hungarian(cost)
    
    #########################################
    # apply the permutation
    #########################################
    costnew = 0.
    X2old = np.copy(X2)
    for (iold, inew) in newind:
        costnew += cost[iold, inew]
        if iold != inew:
            atomiold = atomlist[iold]
            atominew = atomlist[inew]
            X2[atominew,:] = X2old[atomiold,:]


    dist = np.sqrt(costnew)
    return dist, X1.reshape(-1), X2.reshape(-1)

def findBestPermutation( X1, X2, permlist = None ):
    """
    find the permutation of the atoms which minimizes the distance |X1-X2|
    
    Parameters
    ----------
    X1, X2 : 
        the structures to align
    permlist : a list of lists
        A list of lists of atoms which are interchangable.
        e.g. for a 50/50 binary mixture, 
        
            permlist = [range(1,natoms/2), range(natoms/2,natoms)]

    """
    if permlist is None:
        permlist = [range(len(X1)/3)]
    for atomlist in permlist:
        dist, X1, X2 = findBestPermutationList( X1, X2, atomlist )
    dist = np.linalg.norm(X1-X2)
    return dist, X1, X2



import unittest
from testmindist import TestMinDist
class TestMinDistUtils(TestMinDist):
    def setUp(self):
        from pygmin.potentials.ljpshiftfast import LJpshift as BLJ
        from pygmin import defaults
        
        self.natoms = 500
        self.ntypeA = int(self.natoms * .8)
        self.pot = BLJ(self.natoms, self.ntypeA)
        self.permlist = [range(self.ntypeA), range(self.ntypeA, self.natoms)]
        
        self.X1 = np.random.uniform(-1,1,[self.natoms*3])*(float(self.natoms))**(1./3)/2
        ret = defaults.quenchRoutine(self.X1, self.pot.getEnergyGradient, tol=.1)
        self.X1 = ret[0]


#    def testBLJ(self):
#        X1 = np.copy(self.X1)
#        X2 = np.random.uniform(-1,1,[self.natoms*3])*(float(self.natoms))**(1./3)/2
#        
#        #run a quench so the structure is not crazy
#        ret = quench(X2, self.pot.getEnergyGradient)
#        X2 = ret[0]
#
#        self.runtest(X1, X2, minPermDistStochastic)


    def testBLJ_isomer(self):
        """
        test with BLJ potential.  We have two classes of permutable atoms  
        
        test case where X2 is an isomer of X1.
        """
        X1i = np.copy(self.X1)
        X1 = np.copy(self.X1)        
        X2 = np.copy(X1)
        
        #permute X2
        import random, copy
        for atomlist in self.permlist:
            perm = copy.copy(atomlist)
            random.shuffle( perm )
            X2 = permuteArray( X2, perm)

        X2i = np.copy(X2)
        
        #distreturned, X1, X2 = self.runtest(X1, X2)
        distreturned, X1, X2 = self.runtest(X1, X2, findBestPermutation)
        #X1 = X1i
        #X2 = X2i
        #distreturned, X1, X2 = self.runtest(X1, X2, findBestPermutation2)

        
        #it's an isomer, so the distance should be zero
        self.assertTrue( abs(distreturned) < 1e-14, "didn't find isomer: dist = %g" % (distreturned) )

def run_blj():
    from pygmin.potentials.ljpshiftfast import LJpshift as BLJ
    from pygmin import defaults
    
    natoms = 5000
    ntypeA = int(natoms * .8)
    pot = BLJ(natoms, ntypeA)
    permlist = [range(ntypeA), range(ntypeA, natoms)]
    
    X1 = np.random.uniform(-1,1,[natoms*3])*(float(natoms))**(1./3)/2
#    ret = defaults.quenchRoutine(X1, pot.getEnergyGradient, tol=.1)
#    X1 = ret[0]

    X1i = np.copy(X1)
    X1 = np.copy(X1)        
    X2 = np.copy(X1)
    
    #permute X2
    import random, copy
    for atomlist in permlist:
        perm = copy.copy(atomlist)
        random.shuffle( perm )
        X2 = permuteArray( X2, perm)

    X2i = np.copy(X2)
    
    #distreturned, X1, X2 = runtest(X1, X2)
    distreturned, X1, X2 = findBestPermutation(X1, X2, permlist)
    #X1 = X1i
    #X2 = X2i
    #distreturned, X1, X2 = runtest(X1, X2, findBestPermutation2)

    
    #it's an isomer, so the distance should be zero
    if not abs(distreturned) < 1e-14: 
        "didn't find isomer: dist = %g" % (distreturned)

if __name__ == "__main__":
    run_blj()
    
    #unittest.main()
    
    
    

