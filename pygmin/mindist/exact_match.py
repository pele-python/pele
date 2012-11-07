import numpy as np
from mindistutils import findBestPermutation, CoMToOrigin
import rmsfit

class PerfectMatch(object):
    ''' Deterministic check if 2 clusters are a perfect match
    
        Determines quickly if 2 clusters are a perfect match. The algorithm
        does the following:
        
        1) get list of 1st matching candidates in most outer shell
        2) if 1st candidate is already uniquely determined, look 
           in second shell outer for more candidates
        3) loop over all candidate combinations to determine
           orientation and check for match
           
        Notes
        -----
        - The current implementation will fail if the 2 farest atoms are on
        a line
    '''
    
    def __init__(self, accuracy=0.01, permlist=None, align_com = True):
        self.accuracy = accuracy
        self.align_com = align_com
        self.permlist=permlist
        
    def __call__(self, x1, x2):
        natoms = x1.shape[0]
        
        if(self.align_com):
            x1 = CoMToOrigin(x1).reshape([-1,3])
            x2 = CoMToOrigin(x2).reshape([-1,3])
        
        # calculate distance of all atoms
        R1 = np.sqrt(np.sum(x1*x1, axis=1))
        R2 = np.sqrt(np.sum(x1*x1, axis=1))
        
        # at least 2 atoms are needed
        # get atom most outer atom
        R1max, R2max = R1.max(), R2.max()
        
        # do a very quick check
        if np.abs(R1max - R2max) > self.accuracy:
            return False
        
        # get atoms in outer shell
        candidates1_1 = np.arange(natoms)[R1 > R1max - self.accuracy]
        candidates2_1 = np.arange(natoms)[R2 > R2max - self.accuracy]
        
        candidates1_2 = candidates1_1
        # need to check second shell for more candidates?
        if(len(candidates1_2) < 2):
            R1max2 = R1[R1 < R1max].max()
            candidates1_2 = np.arange(natoms)[(R1 > R1max2 - self.accuracy) * (R1 < R1max)]
        
        candidates2_2 = candidates2_1
        # need to check second shell for more candidates?
        if(len(candidates2_2) < 2):
            R2max2 = R2[R2 < R2max].max()
            candidates2_2 = np.arange(natoms)[(R2 > R2max2 - self.accuracy) * (R2 < R2max)]
         
        # now loop over all cobinations
        for i1_1 in candidates1_1:
            for i1_2 in candidates1_2:
                if i1_1 == i1_2:
                    continue
                for i2_1 in candidates2_1:
                    for i2_2 in candidates2_2:
                        if i2_1 == i2_2:
                            continue
                        rot = rmsfit.findrotation_kabsch( \
                                      x1[[i1_1, i1_2]], x2[[i2_1, i2_2]], False)
                        if self.check_match(x1, x2, rot):
                            return True
                        
        return False
    
    def check_match(self, x1, x2, rot):
        x1_trial = np.dot(rot, x1.transpose()).transpose()
        x2_trial = x2.copy()
        dist, x1n, x2n = findBestPermutation(x1_trial, x2_trial)
        max_dist = np.sqrt(np.sum((x1n - x2n)*(x1n - x2n), axis=1)).max()
        if  max_dist < self.accuracy:
            return True
        print max_dist
        return False
                        
    
if __name__ == '__main__':
    natoms = 100
    from pygmin.utils import rotations
    
    for i in xrange(100):
        xx1 = np.random.random(3*natoms)*5
        xx1 = xx1.reshape([-1,3])
        mx = rotations.q2mx(rotations.random_q())
        xx2 = np.dot(mx, xx1.transpose()).transpose()
        xx2 +=2.*(np.random.random(xx2.shape)-0.5)*0.01

        print i,PerfectMatch()(xx1, xx2)
    