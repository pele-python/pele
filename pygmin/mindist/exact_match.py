import numpy as np
from mindistutils import findBestPermutation, CoMToOrigin
import rmsfit


__all__= ["ExactMatchCluster"]

class ExactMatchCluster(object):
    ''' Deterministic check if 2 clusters are a perfect match
    
        Determines quickly if 2 clusters are a perfect match. It uses the 
        atoms which are far away from the center to determine possible
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
    
    def __init__(self, accuracy=0.01, check_inversion=True, permlist=None, align_com = True):
        self.accuracy = accuracy
        self.align_com = align_com
        self.permlist = permlist
        self.check_inversion = check_inversion
        
    def __call__(self, x1, x2):
        if(self.align_com):
            x1 = CoMToOrigin(x1).reshape([-1,3])
            x2 = CoMToOrigin(x2).reshape([-1,3])
        
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
        if np.abs(R1[idx1_1] - R2.max()) > self.accuracy:
            return False
        
        # get indices of atoms in shell of thickness 2*accuracy
        candidates1 = np.arange(len(R2))[ \
             (R2 > R1[idx1_1] - self.accuracy)*(R2 < R1[idx1_1] + self.accuracy)] 
        candidates2 = np.arange(len(R2))[ \
             (R2 > R1[idx1_2] - self.accuracy)*(R2 < R1[idx1_2] + self.accuracy)] 
        
        # now loop over all cobinations
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
                              x1[[idx1_1, idx1_2]], x2[[idx2_1, idx2_2]], False)
                # pass on the match for a closer check
                if self.check_match(x1, x2, rot, 1.0):
                    return True
                
                if self.check_inversion:
                    # get rotation for current atom match candidates
                    rot = rmsfit.findrotation_kabsch( \
                                  x1[[idx1_1, idx1_2]], -x2[[idx2_1, idx2_2]], False)
                    # pass on the match for a closer check
                    if self.check_match(x1, -x2, rot, -1.0):
                        return True
                
                        
        return False
    
    def check_match(self, x1, x2, rot, inverse):
        ''' Make a more detailed comparison if the 2 structures match
        
        Parameters
        ----------
        
        x1: np.array
            coordinates of structure 1
        x2: np.array
            coordinates of structure 2
        rot: np.array, 3x3
            guessed rotation based on reference atoms         
        inverse: double
            -1.0 if do match for inverted coordinates, 1.0 otherwise
            
        returns: boolean
            True or False for match
            
        '''    
        # apply the rotation
        x1_trial = np.dot(rot, x1.transpose()).transpose()
        # make a copy since findBestPermutation will mess up order
        x2_trial = x2.copy()
        # get the best permutation
        dist, x1n, x2n = findBestPermutation(x1_trial.flatten(), x2_trial.flatten())
        x1n = x1n.reshape([-1,3])
        x2n = x2n.reshape([-1,3])
        
        #x1n = x1_trial
        #x2n = x2_trial
        # now find best rotational alignment, this is more reliable than just
        # aligning the 2 reference atoms
        rot2 = rmsfit.findrotation_kabsch(x1n, x2n)
        x1n = np.dot(rot2, x1n.transpose()).transpose()
        
        # use the maximum distance, not rms as cutoff criterion
        max_dist = np.sqrt(np.sum((x1n - x2n)*(x1n - x2n), axis=1)).max()
        if  max_dist < self.accuracy:
            self.best_rotation = np.dot(rot2, rot)
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
        xx2 +=2.*(np.random.random(xx2.shape)-0.5)*0.003
        #xx2 = xx1.copy()
        tmp = xx2[1].copy()
        xx2[1] = xx2[4]
        xx2[4] = tmp
        #dist, x1n, x2n = findBestPermutation(xx1.flatten(), xx2.flatten())
        #print dist
        print i,ExactMatchCluster()(xx1, xx2)
    