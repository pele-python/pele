from mindistutils import getAlignRotation
from permutational_alignment import find_best_permutation
from pygmin.utils import rotations
import numpy as np

class TransformPolicy(object):
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
        ''' translate the coordinates '''
        raise NotImplementedError
    
    def rotate(self, X, mx):
        ''' apply rotation matrix mx for a rotation around the origin'''
        raise NotImplementedError
    
    def can_invert(self):
        ''' returns True or False if an inversion can be performed'''
        raise NotImplementedError
    
    def invert(self, X):
        ''' perform an inversion at the origin '''
        raise NotImplementedError
    
    def permute(self, X, perm):
        ''' returns the permuted coordinates '''
    
class MeasurePolicy(object):
    ''' interface for possible measurements on a set of coordinates
    
    The MeasurePolicy defines an interface which defines how to perform
    certain measures which are essential for minpermdist on a set of
    coordinates. For more motivation of this class see TransformPolicy.
    '''
    
    def get_com(self, X):
        ''' calculate the center of mass '''
        raise NotImplementedError
    
    def get_dist(self, X1, X2):
        ''' calculate the distance between 2 set of coordinates '''
        raise NotImplementedError
    
    def find_permutation(self, X1, X2):
        ''' find the best permutation between 2 sets of coordinates '''
        raise NotImplementedError
    
    def find_rotation(self, X1, X2):
        ''' find the best rotation matrix to bring structure 2 on 1 '''
        raise NotImplementedError

class TransformAtomicCluster(TransformPolicy):
    ''' transformation rules for atomic clusters '''
        
    def translate(self, X, d):
        Xtmp = X.reshape([-1,3])
        Xtmp+=d
        
    def rotate(self, X, mx,):
        Xtmp = X.reshape([-1,3])
        Xtmp = np.dot(mx, Xtmp.transpose()).transpose()
        X[:] = Xtmp.reshape(X.shape)
        
    def permute(self, X, perm):
        return X.reshape(-1,3)[perm].flatten()
        
    def can_invert(self):
        return True
    
    def invert(self, X):
        X[:] = -X
        
class MeasureAtomicCluster(MeasurePolicy):
    ''' measure rules for atomic clusters '''
    
    def __init__(self, permlist=None):
        self.permlist = permlist
    
    def get_com(self, X):
        X = np.reshape(X, [-1,3])
        natoms = len(X[:,0])
        com = X.sum(0) / natoms
        return com

    def get_dist(self, X1, X2):
        return np.linalg.norm(X1-X2)
    
    def find_permutation(self, X1, X2):
        return find_best_permutation(X1, X2, self.permlist)
    
    def find_rotation(self, X1, X2):
        dist, Q2 = getAlignRotation(X1, X2)
        mx = rotations.q2mx(Q2)
        return dist, mx
    
    