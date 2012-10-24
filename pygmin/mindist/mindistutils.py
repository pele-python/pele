import numpy as np
import copy
import pygmin.utils.rotations as rot
import itertools

__all__ = ["alignCoM", "CoMToOrigin", "getAlignRotation", "alignRotation", 
           "findBestPermutation", "findBestPermutationRBMol", "aa2xyz",
           "getDistxyz", "getDistaa", "MinDistWrapper"]

class MinDistWrapper(object):
    """
    wrap a mindist routine into a callable object with the form mindist(X1, X2)
    
    Parameters
    ----------
    mindist : callable
        the mindist routine
    args : 
        extra arguements for mindist
    kwargs : 
        extra keyword arguments for mindist
    """
    def __init__(self, mindist, *args, **kwargs):
        self.mindist = mindist
        self.args = args
        self.kwargs = kwargs
    
    def __call__(self, X1, X2):
        return self.mindist(X1, X2, *self.args, **self.kwargs)

def alignCoM( X1, X2):
    """
    align the center of mass of X2 with that of X1
    """
    natoms = len(X1) / 3
    for i in xrange(3):
        com = np.sum( X1[i::3] - X2[i::3] )
        com /= natoms
        X2[i::3] += com

def CoMToOrigin( X1):
    """
    move the center of mass to the origin
    """
    X1 = np.reshape(X1, [-1,3])
    natoms = len(X1[:,0])
    com = X1.sum(0) / natoms
    X1 -= com
    return X1.reshape(-1)


def getAlignRotation(XA, XB):
    """
    Return the quaternion which aligns XB with XA
    
    Return the quaternion which
    aligns structure XB to be as similar as possible to structure XA.
    To be precise, rotate XB, so as to minimize the distance |XA - XB|.

    Rotations will be done around the origin, not the center of mass

    Rotational alignment follows the prescription of
    Kearsley, Acta Cryst. A, 45, 208-210, 1989
    http://dx.doi.org/10.1107/S0108767388010128
    """
    nsites = len(XA)/3

    #move the center of mass to the origin
    #CoMToOrigin(XA)
    #CoMToOrigin(XB)

    #distcm = np.linalg.norm(XA-XB)
    #print "distcm", distcm

    #########################################
    #Create matrix QMAT
    #########################################

    QMAT = np.zeros([4,4], np.float64)
    for J1 in xrange(nsites):
        J2 = 3*(J1) -1
        XM = XA[J2+1] - XB[J2+1]
        YM = XA[J2+2] - XB[J2+2]
        ZM = XA[J2+3] - XB[J2+3]
        XP = XA[J2+1] + XB[J2+1]
        YP = XA[J2+2] + XB[J2+2]
        ZP = XA[J2+3] + XB[J2+3]
        QMAT[0,0] = QMAT[0,0] + XM**2 + YM**2 + ZM**2
        QMAT[0,1] = QMAT[0,1] - YP*ZM + YM*ZP
        QMAT[0,2] = QMAT[0,2] - XM*ZP + XP*ZM
        QMAT[0,3] = QMAT[0,3] - XP*YM + XM*YP
        QMAT[1,1] = QMAT[1,1] + YP**2 + ZP**2 + XM**2
        QMAT[1,2] = QMAT[1,2] + XM*YM - XP*YP
        QMAT[1,3] = QMAT[1,3] + XM*ZM - XP*ZP
        QMAT[2,2] = QMAT[2,2] + XP**2 + ZP**2 + YM**2
        QMAT[2,3] = QMAT[2,3] + YM*ZM - YP*ZP
        QMAT[3,3] = QMAT[3,3] + XP**2 + YP**2 + ZM**2

    QMAT[1,0] = QMAT[0,1]
    QMAT[2,0] = QMAT[0,2]
    QMAT[2,1] = QMAT[1,2]
    QMAT[3,0] = QMAT[0,3]
    QMAT[3,1] = QMAT[1,3]
    QMAT[3,2] = QMAT[2,3]

    ###########################################
    """
    Find eigenvalues and eigenvectors of QMAT.  The eigenvector corresponding
    to the smallest eigenvalue is the quaternion which rotates XB into best
    alignment with XA.  The smallest eigenvalue is the squared distance between
    the resulting structures.
    """
    ###########################################
    (eigs, vecs) = np.linalg.eig(QMAT)
    #print "eigenvalues", eigs

    imin = np.argmin(eigs)
    eigmin = eigs[imin] #the minimum eigenvector
    Q2 = vecs[:,imin]  #the eigenvector corresponding to the minimum eigenvalue
    if eigmin < 0.:
        if abs(eigmin) < 1e-6:
            eigmin = 0.
        else:
            print 'minDist> WARNING minimum eigenvalue is ',eigmin,' change to absolute value'
            eigmin = -eigmin

    dist = np.sqrt(eigmin) #this is the minimized distance between the two structures
    #print "dist from eigenvalue", dist
    #print "Q2", Q2, "norm", np.linalg.norm(Q2)
    #aa = rot.q2aa( Q2)
    #print "aa ", aa, "norm", np.linalg.norm(aa)

    return dist, Q2


def alignRotation(XA, XB):
    """
    Align structure XB with XA

    Align structure XB to be as similar as possible to structure XA.
    To be precise, rotate XB, so as to minimize the distance |XA - XB|.

    Rotations will be done around the origin, not the center of mass
    """
    nsites = len(XA)/3

    dist, Q2 = getAlignRotation(XA, XB)
    ###################################################################
    # Q2 contains the quaternion which rotates XB to best align with X1.
    # rotate XB according to Q2
    ###################################################################

    rot_mx = rot.q2mx( Q2 )
    #print rot_mx
    for j in range(nsites):
        i = 3*j
        XB[i:i+3] = np.dot( rot_mx, XB[i:i+3] )
    
    dist = np.linalg.norm(XA - XB) #more precise than what it returned

    return dist, XB


def permuteArray(Xold, perm):
    #don't modify Xold
    Xnew = np.copy(Xold)
    permsorted = sorted(perm)
    for (iold, inew) in itertools.izip(permsorted, perm):
        #print iold, "->", inew
        Xnew[inew*3:inew*3+3] = Xold[iold*3:iold*3+3]

    return Xnew



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
    atomlistnp = np.array(atomlist)
    X13 = np.reshape(X1, [-1,3])[atomlistnp,:]
    X23 = np.reshape(X2, [-1,3])[atomlistnp,:]
    cost = (((X13[np.newaxis,:] - X23[:,np.newaxis,:])**2).sum(2))
    #cost = np.sqrt(cost)


    #########################################
    # run the hungarian algorithm
    #########################################
    try:
        #use the hungarian package which is compiled
        import hungarian
        newind1 = hungarian.lap(cost)
        newind = [(i, j) for i,j in enumerate(newind1[0])]
        #print "hungari newind", newind
    except ImportError:
        try:
            #use the munkres package
            #convert cost matrix to a form used by munkres
            from munkres import Munkres
            matrix = cost.tolist()
            m = Munkres()
            newind = m.compute(matrix)
            #print "munkres newind", newind
        except ImportError:
            print "ERROR: findBestPermutation> You must install either the hungarian or the munkres package to use the Hungarian algorithm"
            dist = np.linalg.norm( X1 - X2 )
            return dist, X1, X2

            


    #########################################
    # apply the permutation
    #########################################
    costnew = 0.;
    X2old = np.copy(X2)
    for (iold, inew) in newind:
        costnew    += cost[iold,inew]
        if iold != inew:
            atomiold = atomlist[iold]
            atominew = atomlist[inew]
            #print atomiold, "->", atominew, (iold, inew), "matrix %10.4f, %10.4f" % (matrix[iold][inew], matrix[inew][iold])
            #for i in [iold, inew]:
                #for j in [iold, inew]:
                    #r = np.linalg.norm( X1[i*3:i*3+3] - X2old[j*3:j*3+3] )
                    #print "    %4d %4d %10.4f, %10.4f" % (i, j, r, r**2 ), cost[j, i], matrix[i][j]
            X2[atominew*3:atominew*3+3] = X2old[atomiold*3:atomiold*3+3]

    #costold = sum( [matrix[i][i] for i in range(nsites)] )
    #print "costold    ", costold, np.sqrt(costold)
    #print "costnew    ", costnew, np.sqrt(costnew)

    dist = np.sqrt(costnew)
    return dist, X1, X2

def findBestPermutation( X1, X2, permlist = [] ):
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
    if len(permlist) == 0:
        permlist = [range(len(X1)/3)]
    for atomlist in permlist:
        dist, X1, X2 = findBestPermutationList( X1, X2, atomlist )
    dist = np.linalg.norm(X1-X2)
    return dist, X1, X2

#def findBestPermutation2( X1, X2, permlist = [] ):
#    """
#    find the permutation of the atoms which minimizes the distance |X1-X2|
#    """
#    if len(permlist) == 0:
#        permlist = [range(len(X1)/3)]
#    for atomlist in permlist:
#        dist, X1, X2 = findBestPermutationList2( X1, X2, atomlist )
#    dist = np.linalg.norm(X1-X2)
#    return dist, X1, X2


def molmolMinSymDist(com1, aa1, com2, aa2, mol):
    """
    return the minimal distance between two molecules
    defined by mol and com1, aa1 and mol, com2, aa2
    
    Don't rotate or translate, just apply inner-molecular symmetry operations
    """
    xyz1 = mol.getxyz(com1, aa1)
    mindist = 10000.
    for xyz2, newaa in mol.getSymmetries(com2, aa2):
        dist = np.linalg.norm(xyz1 - xyz2)
        if dist  < mindist:
            mindist = dist
            aamin = newaa.copy()
    return mindist, aamin

def getDistxyz( xyz1, xyz2 ):
    return np.linalg.norm(xyz1 - xyz2)

def getDistaa(coords1, coords2, mysys):
    xyz1 = mysys.getxyz(coords1)
    xyz2 = mysys.getxyz(coords2)
    return getDistxyz(xyz1, xyz2)


def findBestPermutationRBMol_list(coords1, coords2, mol, mollist):
    """
    find the permutation of the molecules which minimizes the distance between the two coordinates
    """    
    nmol = len(coords1) / 3 / 2
    nperm = len(mollist)
    coords2old = coords2.copy()
    #########################################
    # create the cost matrix
    #########################################
    cost = np.zeros( [nperm,nperm], np.float64)
    for i in range(nperm):
        imol = mollist[i]
        com1 = coords1[         imol*3 :          imol*3 + 3]
        aa1  = coords1[3*nmol + imol*3 : 3*nmol + imol*3 + 3]
        for j in range(nperm):
            jmol = mollist[j]
            com2 = coords2[         jmol*3 :          jmol*3 + 3]
            aa2  = coords2[3*nmol + jmol*3 : 3*nmol + jmol*3 + 3]
            cost[j,i], newaa = molmolMinSymDist(com1, aa1, com2, aa2, mol)
    
    #convert cost matrix to a form used by munkres
    matrix = cost.tolist()

    #########################################
    # run the hungarian algorithm
    #########################################
    try:
        #use the hungarian package which is compiled
        import hungarian
        newind1 = hungarian.lap(cost)
        newind = [(i, j) for i,j in enumerate(newind1[0])]
        #print "hungari newind", newind
    except ImportError:
        try:
            #use the munkres package
            #convert cost matrix to a form used by munkres
            from munkres import Munkres
            matrix = cost.tolist()
            m = Munkres()
            newind = m.compute(matrix)
            #print "munkres newind", newind
        except ImportError:
            print "ERROR: findBestPermutation> You must install either the hungarian or the munkres package to use the Hungarian algorithm"
            #raise Exception("ERROR: findBestPermutation> You must install either the hungarian or the munkres package to use the Hungarian algorithm")
            dist = np.linalg.norm( coords1 - coords2 )
            return dist, coords1, coords2


    #########################################
    # apply the permutation
    #########################################
    costnew = 0.;
    coords2 = coords2old.copy()
    for (iold, inew) in newind:
        costnew    += cost[iold,inew]
        if iold != inew:
            moliold = mollist[iold]
            molinew = mollist[inew]
            #print "%4d  ->  %4d" % (moliold, molinew) 
            #change the com coords
            coords2[         molinew*3 : molinew*3+3] = coords2old[         moliold*3 : moliold*3+3]
            #change the aa coords
            coords2[3*nmol + molinew*3 : 3*nmol + molinew*3+3] = coords2old[3*nmol + moliold*3 : 3*nmol + moliold*3+3]

    dist = np.sqrt(costnew)
    return dist, coords1, coords2

def findBestPermutationRBMol(coords1, coords2, mysys, permlist):
    """
    find the best permutation of the molecules.  
    
    then, for each molecule, apply the symmetry operation which minimized the distance
    """
    nmol = len(coords1) / 3 / 2
    if len(permlist) == 0:
        permlist = [range(nmol)]
    for mollist in permlist:
        mol = mysys.molecule_list[ mollist[0] ] #a molecule object corresponding to this permuation
        dist, coords1, coords2 = findBestPermutationRBMol_list( coords1, coords2, mol, mollist )
    dist = getDistaa(coords1, coords2, mysys) 
    #print "in findBest after perm", dist, mysys.getEnergy(coords2)

    #the molecules now have their correct permutation.  
    #For each molecule, apply the symmetry operation which minimized the distance
    for i, mol in enumerate(mysys.molecule_list):
        kcom = 3*i
        kaa  = 3*nmol + 3*i
        com1 = coords1[kcom : kcom + 3]
        aa1  = coords1[kaa  : kaa  + 3]
        com2 = coords2[kcom : kcom + 3]
        aa2  = coords2[kaa  : kaa  + 3]
        d, aanew = molmolMinSymDist(com1, aa1, com2, aa2, mol)
        coords2[kaa : kaa + 3] = aanew
    dist = getDistaa(coords1, coords2, mysys) 
    #print "after inner-mol sym", dist
    return dist, coords1, coords2


    
  


def aa2xyz(XB, AA):
    """
    Rotate XB according to angle axis AA
    """
    nsites = len(XB)/3
    XBnew = np.copy(XB)
    rot_mx = rot.aa2mx( AA )
    for j in range(nsites):
        i = 3*j
        XBnew[i:i+3] = np.dot( rot_mx, XBnew[i:i+3] )
    return XBnew


import unittest
from testmindist import TestMinDist
class TestMinDistUtils(TestMinDist):
    def setUp(self):
        from pygmin.potentials.ljpshift import LJpshift as BLJ
        from pygmin import defaults
        
        self.natoms = 10
        self.ntypeA = int(self.natoms * .8)
        self.pot = BLJ(self.natoms, self.ntypeA)
        self.permlist = [range(self.ntypeA), range(self.ntypeA, self.natoms)]
        
        self.X1 = np.random.uniform(-1,1,[self.natoms*3])*(float(self.natoms))**(1./3)/2
        ret = defaults.quenchRoutine(self.X1, self.pot.getEnergyGradient, **defaults.quenchParams)
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

if __name__ == "__main__":
    unittest.main()
