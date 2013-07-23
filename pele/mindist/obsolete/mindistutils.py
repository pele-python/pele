import numpy as np
import copy
import pele.utils.rotations as rot
import itertools

__all__ = ["alignCoM", "CoMToOrigin", "getAlignRotation", "alignRotation", 
           "findBestPermutationRBMol", "aa2xyz",
           "getDistxyz", "getDistaa"]

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
    if permlist is None:
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


