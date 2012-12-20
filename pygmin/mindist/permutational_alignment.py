import numpy as np
import itertools

__all__ = ["findBestPermutation"] 

have_minperm = False
have_hungarian = False
have_munkres = False
try:
    import minperm
    have_minperm = True    
except ImportError:
    pass
try:
    import hungarian
    have_hungarian = True
except ImportError:
    pass
try:    
    import munkres
    have_munkres = True
except ImportError:
    pass


_findBestPermutationList = None
if have_minperm:
    def _findBestPermutationList(*args, **kwargs):
        return findBestPermutationListOPTIM(*args, **kwargs)

elif have_hungarian:
    def _findBestPermutationList(*args, **kwargs):
        return findBestPermutationListHungarian(*args, **kwargs)
elif have_munkres:
    def _findBestPermutationList(*args, **kwargs):
        return findBestPermutationListMunkres(*args, **kwargs)
else:
    raise BaseException("No Hungarian algorithm implementation found!"
                        "Please compile minperm.f90 or install the hungarian or the munkres package")





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

def findBestPermutationListMunkres( X1, X2, atomlist = None ):
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

    #########################################
    # create the cost matrix
    # cost[j,i] = (X1(i,:) - X2(j,:))**2
    #########################################
    X1 = X1.reshape([-1,3])
    X2 = X2.reshape([-1,3])
    cost = makeCostMatrix(X1, X2, atomlist)
    #cost = np.sqrt(cost)

    #########################################
    # run the munkres algorithm
    #########################################
    matrix = cost.tolist()
    m = munkres.Munkres()
    newind = m.compute(matrix)    
    
    #########################################
    # apply the permutation
    #########################################
    costnew = 0.
    X2new = np.copy(X2)
    for (iold, inew) in newind:
        costnew += cost[iold, inew]
        if iold != inew:
            atomiold = atomlist[iold]
            atominew = atomlist[inew]
            X2new[atominew,:] = X2[atomiold,:]
        
    X1 = X1.reshape(-1)
    X2new = X2new.reshape(-1)
#    dist = np.linalg.norm(X1-X2new)
    dist = np.sqrt(costnew)
    return dist, X1.reshape(-1), X2new.reshape(-1)

def findBestPermutationListHungarian( X1, X2, atomlist = None ):
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
    atomlist = np.array(atomlist)


    #########################################
    # create the cost matrix
    # cost[j,i] = (X1(i,:) - X2(j,:))**2
    #########################################
    X1 = X1.reshape([-1,3])
    X2 = X2.reshape([-1,3])
    cost = makeCostMatrix(X1, X2, atomlist)
    #cost = np.sqrt(cost)

    #########################################
    # run the hungarian algorithm
    #########################################
    newind1 = hungarian.lap(cost)
    perm = newind1[1]

    #note: the hungarian algorithm changes
    #the cost matrix.  I'm not sure why, and it may be a bug, 
    #but the indices it returns are still correct
#    if not np.all(cost >= 0):
#        m = np.max(np.abs(cost-costsave))
#        print "after hungarian cost greater than zero:, %g" % m
    
    
    #########################################
    # apply the permutation
    #########################################
    newperm = np.array(atomlist[perm])
    X2new = np.copy(X2)
    X2new[atomlist,:] = X2[newperm,:]

        
    X1 = X1.reshape(-1)
    X2new = X2new.reshape(-1)
    dist = np.linalg.norm(X1-X2new)
    return dist, X1, X2new

def findBestPermutationListOPTIM(X1, X2, atomlist, boxl=None):
    """
    use OPTIM's minperm() routine to calculate the optimum permutation
    """    
    #deal with periodic boundary conditions
    periodic = boxl is not None
    if not periodic:
        #it must have a value for passing to fortran 
        boxl = 1.
    sx = sy = sz = boxl
    
    X1 = X1.reshape([-1,3])
    X2 = X2.reshape([-1,3])
    atomlist = np.array(atomlist)
    X13 = X1[atomlist,:].reshape(-1)
    X23 = X2[atomlist,:].reshape(-1)
    
    #run the minperm algorithm
    perm, dist, worstdist, worstradius = minperm.minperm(X13, X23, sx, sy, sz, periodic)
    perm -= 1 #fortran indexing

    #note, dist returned by minperm comes will only be accurate to 6 decimal places at best.
    #if we want a more accurate distance we should calculate it from the coordinates

    # apply the permutation
    newperm = np.array(atomlist[perm])
    X2new = np.copy(X2)
    X2new[atomlist,:] = X2[newperm,:]

    dist = np.sqrt(dist)
    return dist, X1.reshape(-1), X2new.reshape(-1)

    

def findBestPermutation( X1, X2, permlist = None, user_algorithm=None):
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
    user_algoriithm : None or callable
        you can optionally pass which algorithm to use.
    
    Returns
    -------
    dist : float
        the minimum distance
    X1new : 
        should be the same as X1
    X2new :
        X2 in best alignment with X1
    
    Notes
    -----
    This for each list of interchangeable atoms in permlist the permutation
    which minimizes the distance between the two structures is found.  This minimimization
    is done by mapping the problem onto the linear assignment problem which can then be solved
    using graph theoretic techniques.  
    
    http://en.wikipedia.org/wiki/Linear_assignment_problem
    http://en.wikipedia.org/wiki/Hungarian_algorithm

    there are several packages in pypi which solve the linear assignment problem
    
    hungarian : c++ code wrapped in python.  scales roughly like natoms**2.5
    
    munkres : completely in python. scales roughly like natoms**3.  very slow for natoms > 10
    
    in addition we have wrapped the OPTIM version for use in pygmin.  It uses the sparse 
    version of the Jonker-Volgenant algorithm.  Furthermore the cost matrix calculated in 
    a compiled language for an additional speed boost. It scales roughtly like natoms**2

    """
    if permlist is None:
        permlist = [range(len(X1)/3)]
    for atomlist in permlist:
        if user_algorithm is None:
            dist, X1, X2 = _findBestPermutationList( X1, X2, atomlist )
        else:
            dist, X1, X2 = user_algorithm( X1, X2, atomlist )
    dist = np.linalg.norm(X1-X2)
    return dist, X1, X2


#####################################################################
# only testing stuff below here
#####################################################################


import unittest
from testmindist import TestMinDist
class TestPermLJ(TestMinDist):
    """
    test permutational optimization algorithms with the LJ potential
    """
    def setUp(self):
        from pygmin.potentials import LJ
        from pygmin import defaults
        from pygmin.mindist import CoMToOrigin
        
        self.natoms = 15
        self.pot = LJ(self.natoms)
        self.permlist = [range(self.natoms)]
        
        self.X1 = np.random.uniform(-1,1,[self.natoms*3])*(float(self.natoms))**(1./3)/2
        ret = defaults.quenchRoutine(self.X1, self.pot.getEnergyGradient, tol=.1)
        self.X1 = ret[0]
        self.X1 = CoMToOrigin(self.X1)

    def testLJ(self):
        """basic test to make sure everythings working right"""
        import pygmin.defaults as defaults
        X1 = np.copy(self.X1)
        X2 = np.random.uniform(-1,1,[self.natoms*3])*(float(self.natoms))**(1./3)/2
        
        #run a quench so the structure is not crazy
        ret = defaults.quenchRoutine(X2, self.pot.getEnergyGradient)
        X2 = ret[0]

        self.runtest(X1, X2, findBestPermutation)

    def testLJ_OPTIM(self):
        """test findBestPermutationListOPTIM"""
        import pygmin.defaults as defaults
        X1 = np.copy(self.X1)
        X2 = np.random.uniform(-1,1,[self.natoms*3])*(float(self.natoms))**(1./3)/2
        
        #run a quench so the structure is not crazy
        ret = defaults.quenchRoutine(X2, self.pot.getEnergyGradient)
        X2 = ret[0]

        self.runtest(X1, X2, findBestPermutationListOPTIM, atomlist=self.permlist[0])

    def testLJ_munkres(self):
        """test findBestPermutationListOPTIM"""
        import pygmin.defaults as defaults
        X1 = np.copy(self.X1)
        X2 = np.random.uniform(-1,1,[self.natoms*3])*(float(self.natoms))**(1./3)/2
        
        #run a quench so the structure is not crazy
        ret = defaults.quenchRoutine(X2, self.pot.getEnergyGradient)
        X2 = ret[0]

        self.runtest(X1, X2, findBestPermutationListMunkres, atomlist=self.permlist[0])

    def testLJ_hungarian(self):
        """test findBestPermutationListOPTIM"""
        import pygmin.defaults as defaults
        X1 = np.copy(self.X1)
        X2 = np.random.uniform(-1,1,[self.natoms*3])*(float(self.natoms))**(1./3)/2
        
        #run a quench so the structure is not crazy
        ret = defaults.quenchRoutine(X2, self.pot.getEnergyGradient)
        X2 = ret[0]

        self.runtest(X1, X2, findBestPermutationListHungarian, atomlist=self.permlist[0])

    def test_multiple(self):
        """test hungarian, munkres, and OPTIM algorithms agains each other"""
        import pygmin.defaults as defaults
        from pygmin.mindist import CoMToOrigin
        X1 = np.copy(self.X1)
        X2 = np.random.uniform(-1,1,[self.natoms*3])*(float(self.natoms))**(1./3)/2
        X2 = CoMToOrigin(X2)
        X1 = CoMToOrigin(X1)
        X2i = X2.copy()
        
        d1, X11, X21 = findBestPermutationListHungarian(X1, X2, self.permlist[0])
        d1calc = np.linalg.norm(X11-X21)
        
        X2 = X2i.copy()
        d2, X12, X22 = findBestPermutationListOPTIM(X1, X2, self.permlist[0])
        d2calc = np.linalg.norm(X12-X22)

        X2 = X2i.copy()
        d3, X13, X23 = findBestPermutationListMunkres(X1, X2, self.permlist[0])
        d3calc = np.linalg.norm(X13-X23)

        
        self.assertAlmostEqual(d1, d2, 5)
        self.assertAlmostEqual(d1, d1calc, 5)
        self.assertAlmostEqual(d2, d2calc, 5)
        
        self.assertAlmostEqual(d1, d3, 5)
        self.assertAlmostEqual(d3, d3calc, 5)
        
        


import unittest
from testmindist import TestMinDist
class TestMinDistUtils(TestMinDist):
    def setUp(self):
        from pygmin.potentials.ljpshiftfast import LJpshift as BLJ
        from pygmin import defaults
        
        self.natoms = 15
        self.ntypeA = int(self.natoms * .8)
        self.pot = BLJ(self.natoms, self.ntypeA)
        self.permlist = [range(self.ntypeA), range(self.ntypeA, self.natoms)]
        
        self.X1 = np.random.uniform(-1,1,[self.natoms*3])*(float(self.natoms))**(1./3)/2
        ret = defaults.quenchRoutine(self.X1, self.pot.getEnergyGradient, tol=.1)
        self.X1 = ret[0]


    def testBLJ(self):
        import pygmin.defaults as defaults
        X1 = np.copy(self.X1)
        X2 = np.random.uniform(-1,1,[self.natoms*3])*(float(self.natoms))**(1./3)/2
        
        #run a quench so the structure is not crazy
        ret = defaults.quenchRoutine(X2, self.pot.getEnergyGradient)
        X2 = ret[0]

        self.runtest(X1, X2, findBestPermutation)


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
        
        #it's an isomer, so the distance should be zero
        self.assertAlmostEqual(distreturned, 0., 14, "didn't find isomer: dist = %g" % (distreturned) )

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
    
    unittest.main()
    
    
    

