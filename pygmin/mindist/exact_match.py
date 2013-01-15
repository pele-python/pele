import numpy as np
from permutational_alignment import find_best_permutation
from _minpermdist_policies import TransformAtomicCluster, MeasureAtomicCluster
import rmsfit

__all__= ["StandardClusterAlignment", "ExactMatchCluster"]

class StandardClusterAlignment(object):
    '''
    class to iterate over standard alignments for atomic clusters
    
    Quickly determines possible alignments of clusters which exactly match.
    It uses atoms which are far away from the center to determine possible
    rotations. The algorithm does the following:
    
    1) Get 2 reference atoms from structure 1 which are farest away from center
       and are not linear
    2) Determine candidates from structure 2 which are in same shell
       as reference atoms from structure 1 (+- accuracy)
    3) loop over all candidate combinations to determine
       orientation and check for match. Skip directly if angle of candidates
       does not match angle of reference atoms in structure 1.
       
    Parameters
    ----------
    coords1 : np.array
        first coordinates
    coords2 : np.array
        second coordinates
    accuracy : float
        accuracy of shell for atom candidates in standard alignment
    can_invert : boolean
        is an inversion possible?
        
    Examples
    --------
    
    >> for rot, invert in StandardClusterAlignment(X1, X2):
    >>     print "possible rotation:",rot,"inversion:",invert
    
    '''
    def __init__(self, coords1, coords2, accuracy = 0.01, can_invert=True):
        x1 = coords1.reshape([-1,3]).copy()
        x2 = coords2.reshape([-1,3]).copy()
        
        self.accuracy = accuracy
        self.can_invert = can_invert
        
        # calculate distance of all atoms
        R1 = np.sqrt(np.sum(x1*x1, axis=1))
        R2 = np.sqrt(np.sum(x2*x2, axis=1))
        
        # at least 2 atoms are needed
        # get atom most outer atom
        
        # get 1. reference atom in configuration 1
        # use the atom with biggest distance to com
        idx_sorted = R1.argsort()
        idx1_1 = idx_sorted[-1]
        
        # find second atom which is not in a line
        for idx1_2 in reversed(idx_sorted[0:-1]):
            # stop if angle is larger than threshold
            cos_theta1 = np.dot(x1[idx1_1], x1[idx1_2]) / \
                (np.linalg.norm(x1[idx1_1])*np.linalg.norm(x1[idx1_2])) 
            if cos_theta1 < 0.9:
                break
            
        # do a very quick check if most distant atom from
        # center are within accuracy
        if np.abs(R1[idx1_1] - R2.max()) > accuracy:
            candidates1 = []
            candidates2 = []
        else:
            # get indices of atoms in shell of thickness 2*accuracy
            candidates1 = np.arange(len(R2))[ \
                 (R2 > R1[idx1_1] - accuracy)*(R2 < R1[idx1_1] + accuracy)] 
            candidates2 = np.arange(len(R2))[ \
                 (R2 > R1[idx1_2] - accuracy)*(R2 < R1[idx1_2] + accuracy)] 
        
        self.x1 = x1
        self.x2 = x2
        self.idx1_1 = idx1_1
        self.idx1_2 = idx1_2
        self.idx2_1 = None
        self.idx2_2 = None
        self.invert = False
        
        self.candidates2 = candidates2
        
        self.iter1 = iter(candidates1)
        self.iter2 = iter(self.candidates2)
                
    def __iter__(self):
        return self
    
    def next(self):
        # obtain first index for first call
        if self.idx2_1 is None:
            self.idx2_1 = self.iter1.next()
            
        # toggle inversion if inversion is possible
        if self.can_invert and self.invert == False and self.idx2_2 is not None:
            self.invert = True
        else:
            # determine next pair of indices
            self.invert = False
            # try to increment 2nd iterator
            try: 
                self.idx2_2 = self.iter2.next()
            except StopIteration:
                # end of list, start over again
                self.iter2 = iter(self.candidates2)
                # and increment iter1
                self.idx2_1 = self.iter1.next()
                self.idx2_2 = None
                return self.next()
            
        if self.idx2_1 == self.idx2_2:
            return self.next()
        
        x1 = self.x1
        x2 = self.x2
        idx1_1 = self.idx1_1
        idx1_2 = self.idx1_2
        idx2_1 = self.idx2_1
        idx2_2 = self.idx2_2
        
        assert idx1_1 is not None
        assert idx1_2 is not None
        assert idx2_1 is not None
        assert idx2_2 is not None
        
        # we can immediately trash the match if angle does not match
        try:
            cos_theta2 = np.dot(x2[idx2_1], x2[idx2_2]) / \
                (np.linalg.norm(x2[idx2_1])*np.linalg.norm(x2[idx2_2]))
        except ValueError:
            raise
        if(np.abs(cos_theta2 - cos_theta2) > 0.5):
            return self.next()

        mul = 1.0
        if(self.invert):
            mul=-1.0

        # get rotation for current atom match candidates
        dist, rot = rmsfit.findrotation( \
                      x1[[idx1_1, idx1_2]], mul*x2[[idx2_1, idx2_2]], align_com=False)
                
        return rot, self.invert
    
class ExactMatchCluster(object):
    ''' Deterministic check if 2 clusters are a perfect match
    
        Determines quickly if 2 clusters are a perfect match. It uses
        check_standard_alignment_cluster to get possible orientations.
        
        
        
        Parameters
        ----------
        
        accuracy: float, optional
            maximum deviation of atoms to still consider cluster as a match
            
        check_inversion: boolean, optional
            check for inversion symmetry, default is True
            
        permlist: iteratable, optional
            list of allowed permutations. Default is None which means all
            particles can be permuted
            
        align_com: boolean, optional
            Flag if com should be removed before comparison
            
        Examples
        --------
        
        >>> x1 = np.random.random(3*natoms)
        >>> x2 = x1 + 1e-4*np.random.random(x1.shape)
        >>> matches = ExactClusterMatch(accuracy=1e-3)
        >>> if match(x1, x2):
        >>>     print "the two structures are identical
          
    '''
    
    def __init__(self, tol = 0.01, accuracy=0.01, transform=TransformAtomicCluster(), measure=MeasureAtomicCluster()):
        self.accuracy = accuracy
        self.tol = tol
        self.transform = transform
        self.measure = measure
        
        
    def __call__(self, coords1, coords2):
        com1 = self.measure.get_com(coords1)
        x1 = coords1.copy()
        self.transform.translate(x1, -com1)
        
        com2 = self.measure.get_com(coords2)
        x2 = coords2.copy()
        self.transform.translate(x2, -com2)
        
        for rot, invert in StandardClusterAlignment(x1, x2, accuracy = self.accuracy,
                                   can_invert=self.transform.can_invert()):
            if self.check_match(x1, x2, rot, invert):
                return True
        return False
                        
    def check_match(self, x1, x2, rot, invert):
        ''' Make a more detailed comparison if the 2 structures match
        
        Parameters
        ----------
        
        rot: np.array, 3x3
            guessed rotation based on reference atoms         
        invert: boolean
            True do match for inverted coordinates
                        
        returns: boolean
            True or False for match
            
        '''    
        # apply the rotation
        x2_trial = x2.copy()
        if(invert):
            self.transform.invert(x2_trial)
        self.transform.rotate(x2_trial, rot)

        
        # get the best permutation
        dist, perm = self.measure.find_permutation(x1, x2_trial)
        x2_trial = self.transform.permute(x2_trial, perm)
       
        # now find best rotational alignment, this is more reliable than just
        # aligning the 2 reference atoms
        dist, rot2 = self.measure.find_rotation(x1, x2_trial)
        self.transform.rotate(x2_trial, rot2)
        # use the maximum distance, not rms as cutoff criterion
        
        if  self.measure.get_dist(x1, x2_trial) < self.tol:
            return True
        return False
                        
    
if __name__ == '__main__':
    natoms = 35
    from pygmin.utils import rotations
    
    for i in xrange(100):
        xx1 = np.random.random(3*natoms)*5
        xx1 = xx1.reshape([-1,3])
        mx = rotations.q2mx(rotations.random_q())
        xx2 = -np.dot(mx, xx1.transpose()).transpose()
        xx2 +=2.*(np.random.random(xx2.shape)-0.5)*0.001
        #xx2 = xx1.copy()
        tmp = xx2[1].copy()
        xx2[1] = xx2[4]
        xx2[4] = tmp
        #dist, x1n, x2n = findBestPermutation(xx1.flatten(), xx2.flatten())
        #print dist
        print i,ExactMatchCluster()(xx1.flatten(), xx2.flatten())
    
