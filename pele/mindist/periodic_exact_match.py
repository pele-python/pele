import numpy as np
import copy

from _minpermdist_policies import MeasurePolicy, TransformPolicy
from pele.mindist.permutational_alignment import find_best_permutation
from pele.utils.rbtools import CoordsAdapter

class MeasurePeriodic(MeasurePolicy):
    ''' interface for possible measurements on a set of coordinates with periodic boundary conditions
    
    Notes
    -----
    this is only implemented for a rectangular box
    
    Parameters
    ----------
    box_lengths : array
        vector defining the box
    permlist : list of lists
        list of lists of identical atoms
    '''
    def __init__(self, box_lengths, permlist=None):
        self.boxlengths = np.array(box_lengths)
        self.iboxlengths = 1. / self.boxlengths
        self.permlist = permlist

    def get_dist(self, X1, X2):
        ''' calculate the distance between 2 set of coordinates '''
    
        dx = X2 - X1
        dx = dx.reshape(-1,len(self.boxlengths))
        dx -= np.round(dx * self.iboxlengths) * self.boxlengths
        return np.linalg.norm(dx.flatten())
    
    def find_permutation(self, X1, X2):
        return find_best_permutation(X1, X2, self.permlist, box_lengths=self.boxlengths)
    
    def get_com(self, X):
        raise NotImplementedError("Center of mass not defined for periodic systems")   
    

# sn402: do we care about atomic positions, or just centres of mass?
# If the latter, ignore this class altogether
class MeasurePeriodicRigid(MeasurePeriodic):
    def __init__(self, box_lengths, topology, permlist=None):
        self.boxlengths = np.array(box_lengths)
        self.iboxlengths = 1. / self.boxlengths
        self.permlist = permlist
        self.topology = topology
        
    def get_dist(self, X1, X2):                         
        x1 = X1.copy()
        x2 = X2.copy()

        atom1 = self.topology.to_atomistic(x1)
        atom2 = self.topology.to_atomistic(x2)
        
        dx = atom2 - atom1
        dx = dx.reshape(-1,len(self.boxlengths))
        dx -= np.round(dx * self.iboxlengths) * self.boxlengths
        return np.linalg.norm(dx.flatten())
    
    
class TransformPeriodic(TransformPolicy):
    ''' interface for possible transformations on a set of coordinates
    
    The transform policy tells minpermdist how to perform transformations, 
    i.e. a translation, rotation and inversion on a specific set of
    coordinates. This class is necessary since in general a coordinate array
    does not carry any information  on the type of coordinate, e.g. if it's a
    site coordinate, atom coordinate or angle axis vector.
    
    All transformation act in place, that means they change the current
    coordinates and do not make a copy.
    
    '''
    
    def translate(self, X, d):
        Xtmp = X.reshape([-1,3])
        Xtmp += d
    
    def can_invert(self):
        ''' returns True or False if an inversion can be performed'''
        return False
    
    def permute(self, X, perm):
        return X.reshape(-1,3)[perm].flatten()
    
    
    # sn402: new class
class TransformPeriodicRigid(TransformPeriodic):
    # sn402: changed this for rigid body system.
    def translate(self,X,d):
        ca = CoordsAdapter(coords=X)
        if(ca.nrigid > 0):
            ca.posRigid += d

        if(ca.natoms > 0):
            ca.posAtom += d



class ExactMatchPeriodic(object):
    """Deterministic check if 2 structures are a perfect match
    
    assume there are permutable atoms, because if there are not then the problem is trivial.
    if even one atom is not permutable, then it is greatly simplified.  if there are two atoms 
    not permutable it becomes trivial.
    
    
    we have two structures, A and B.  
    
    1. Choose an atom iA in structure A.
    
    #. Choose an atom iB in structure B
    
    #. align the structures based on the assumption iA == iB.
    
    #. if the structures are not the same, repeat for all atoms iB in B
    """
    def __init__(self, measure, accuracy=.01):
        self.transform = TransformPeriodic()
        self.measure = measure
        self.permlist = measure.permlist
        self.accuracy = accuracy
    
    def __call__(self, x1, x2):
        x1 = x1.reshape(-1,3)
        x2 = x2.reshape(-1,3)
        x2_init = x2.copy()
        x2 = x2.copy()
        permlist = self.permlist
        
        # get the shortest atomlist from permlist
        if permlist is None:
            atomlist = [range(len(x1.shape[0]))]
        elif len(permlist) == 0:
            # no permutable atoms
            atomlist = [0]
        else:
            atomlist = sorted(permlist, key=lambda a: len(a))[0]

        iA = atomlist[0]
        for iB in atomlist:
            # overlay structures with atom iA == atom iB and check for exact match 
            x2 = x2_init.copy()
            are_match = self.check_match(x1, x2, iA, iB)
            if are_match:
                return True
        
        return False
            
    def check_match(self, x1, x2, iA, iB):
        """overlay structures with atom iA == atom iB and check for exact match""" 
        self.transform.translate(x2, x1[iA,:] - x2[iB,:])
        dist, perm = self.measure.find_permutation(x1, x2)
        x2 = self.transform.permute(x2, perm)
        x2 = x2.reshape(-1,3)
        dist = self.measure.get_dist(x1, x2)
        if dist <= self.accuracy:
            return True
        
        
        
def randomly_permute(x, permlist):
    import random
    x = x.reshape(-1,3)
    xnew = x.copy()
    for atomlist in permlist:
        permutation = copy.copy(atomlist)
        random.shuffle(permutation)
        xnew[atomlist,:] = x[permutation,:]
    return xnew.flatten()

#
# testing only beyond this point
#

        
if __name__ == "__main__":
    from pele.systems import LJCluster
    natoms = 20
    rho = .5
    boxl = (float(natoms) / rho)**(1./3)
    boxlengths = np.ones(3) * boxl
    
    permlist = [range(natoms)]
    measure = MeasurePeriodic(boxlengths, permlist)
    
    system = LJCluster(natoms)

    x1 = system.get_random_configuration()
    x2 = x1.copy()
    x2 = randomly_permute(x2, permlist)
    
    exact_match = ExactMatchPeriodic(measure)
    em = exact_match(x1, x2)
    print em
    
