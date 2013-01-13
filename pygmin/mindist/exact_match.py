import numpy as np
from permutational_alignment import find_best_permutation
from _minpermdist_policies import TransformAtomicCluster, MeasureAtomicCluster
import rmsfit


__all__= ["ExactMatchCluster"]

def check_standard_alignment_cluster(coords1, coords2, check_match, accuracy = 0.01, check_inversion=True):
        '''
        eumerates standard alignments for atomic clusters
        
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
        '''
        x1 = coords1.reshape([-1,3])
        x2 = coords2.reshape([-1,3])
        
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
            return False
        
        # get indices of atoms in shell of thickness 2*accuracy
        candidates1 = np.arange(len(R2))[ \
             (R2 > R1[idx1_1] - accuracy)*(R2 < R1[idx1_1] + accuracy)] 
        candidates2 = np.arange(len(R2))[ \
             (R2 > R1[idx1_2] - accuracy)*(R2 < R1[idx1_2] + accuracy)] 
        
        # now loop over all combinations
        for idx2_1 in candidates1:
            for idx2_2 in candidates2:
                if idx2_1 == idx2_2:
                    continue
                
                # we can immediately trash the match if angle does not match
                cos_theta2 = np.dot(x2[idx2_1], x2[idx2_2]) / \
                    (np.linalg.norm(x2[idx2_1])*np.linalg.norm(x2[idx2_2]))
                if(np.abs(cos_theta2 - cos_theta2) > 0.5):
                    continue

                # get rotation for current atom match candidates
                rot = rmsfit.findrotation_kabsch( \
                              x1[[idx1_1, idx1_2]], x2[[idx2_1, idx2_2]], align_com=False)
                # pass on the match for a closer check
                if check_match(rot.transpose(), False):
                    return True
                
                if check_inversion:
                    # get rotation for current atom match candidates
                    rot = rmsfit.findrotation_kabsch( \
                                  x1[[idx1_1, idx1_2]], -x2[[idx2_1, idx2_2]], align_com=False)
                    # pass on the match for a closer check
                    if check_match(rot.transpose(), True):
                        return True                        
        return False

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
        self.x1 = coords1.copy()
        self.transform.translate(self.x1, -com1)
        
        com2 = self.measure.get_com(coords2)
        self.x2 = coords2.copy()
        self.transform.translate(self.x2, -com2)
        
        return check_standard_alignment_cluster(self.x1, self.x2, self.check_match,
                                   accuracy = self.accuracy,
                                   check_inversion=self.transform.can_invert())
                        
    def check_match(self, rot, invert):
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
        x1 = self.x1
        x2_trial = self.x2.copy()
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
        xx2 = np.dot(mx, xx1.transpose()).transpose()
        xx2 +=2.*(np.random.random(xx2.shape)-0.5)*0.001
        #xx2 = xx1.copy()
        tmp = xx2[1].copy()
        xx2[1] = xx2[4]
        xx2[4] = tmp
        #dist, x1n, x2n = findBestPermutation(xx1.flatten(), xx2.flatten())
        #print dist
        print i,ExactMatchCluster()(xx1.flatten(), xx2.flatten())
    
