import numpy as np
import copy

from _minpermdist_policies import MeasurePolicy, TransformPolicy
from pygmin.mindist import find_best_permutation


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
        Xtmp+=d
    
    def can_invert(self):
        ''' returns True or False if an inversion can be performed'''
        return False
    
    def permute(self, X, perm):
        return X.reshape(-1,3)[perm].flatten()




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
        
        #get the shortest atomlist from permlist
        if permlist is None:
            atomlist = [range(len(x1.shape[0]))]
        elif len(permlist) == 0:
            #no permutable atoms
            atomlist = [0]
        else:
            atomlist = sorted(permlist, key=lambda a: len(a))[0]

        iA = atomlist[0]
        for iB in atomlist:
            # overlay structures with atom iA == atom iB and check for exact match 
            x2 = x2_init.copy()
#            print iA, iB
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

import unittest
class TestExactMatchPeriodicLJ(unittest.TestCase):
    def setUp(self):
        self.natoms = 100
        rho = .5
        boxl = (float(self.natoms) / rho)**(1./3)
        boxlengths = np.ones(3) * boxl + np.random.rand(3)*.1
        
        self.permlist = self.get_permlist()        
        self.measure = MeasurePeriodic(boxlengths, self.permlist)
        
        self.x1 = self.get_random_configuration()
        self.x2diff = self.get_random_configuration()
        self.x2same = self.randomly_permute(self.x1.copy())
        
        self.exact_match = ExactMatchPeriodic(self.measure, accuracy=1e-5)
    
    def get_permlist(self):
        return [range(self.natoms)]
    
    def randomly_permute(self, x):
        import random
        x = x.reshape(-1,3)
        xnew = x.copy()
        for atomlist in self.permlist:
            permutation = copy.copy(atomlist)
            random.shuffle(permutation)
            xnew[atomlist,:] = x[permutation,:]
        return xnew.flatten()

    
    def get_random_configuration(self):
        return np.random.uniform(-1,1,self.natoms*3) * (float(self.natoms))**(1./3) * 1.5
    
    def test_exact_match(self):
        self.assertTrue(self.exact_match(self.x1, self.x2same))
        
    def test_no_exact_match(self):
        self.assertFalse(self.exact_match(self.x1, self.x2diff))

    def test_exact_match_periodic(self):
        self.x2same[:3] += self.measure.boxlengths  
        self.assertTrue(self.exact_match(self.x1, self.x2same))


class TestExactMatchPeriodicBLJ(TestExactMatchPeriodicLJ):        
    def get_permlist(self):
        self.ntypeA = int(self.natoms *.8)
        return [range(self.ntypeA), range(self.ntypeA, self.natoms)]
    
if __name__ == '__main__':
    unittest.main()
        
if __name__ == "__main__":
    from pygmin.systems import LJCluster
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
    em = exact_match.are_exact(x1, x2)
    print em
    
    unittest.main()
    
    
    
    
    
