from __future__ import print_function
import numpy as np
from pele.mindist.periodic_exact_match import TransformPeriodic
from inspect import stack

class MinDistBulk(object):
    """ Obtain the best alignment between two configurations of a periodic system"""
    def __init__(self, boxvec, measure, transform=TransformPeriodic(), niter=10, verbose=False, tol=0.01, 
                 accuracy=0.01):        
        self.niter = niter       
        self.verbose = verbose
        self.measure = measure
        self.transform=transform
        self.accuracy = accuracy
        self.tol = tol
        self.boxvec = boxvec
                          
    def align_fragments(self, x1, x2): 
        """
        Obtain the best alignment between two configurations of a periodic system
        
        Parameters
        ----------
        coords1, coords2 : np.array 
            the structures to align.  X2 will be aligned with X1
            Both structures are arrays of cartesian coordinates

        Returns
        -------
        a triple of (dist, coords1, coords2). coords1 are the unchanged coords1
        and coords2 are brought in best alignment with coords2
        """

        if self.verbose:
            print("Transform:")
            print(self.transform)
            print("Measure.topology:")
            print(self.measure.topology)
            print("Called by", stack())
        
        # we don't want to change the given coordinates
        x1 = np.copy(x1).reshape(-1, self.boxvec.size)
        x2 = np.copy(x2).reshape(-1, self.boxvec.size)

        dist, dx = self.measure.get_dist(x2, x1, with_vector=True)
        dx = dx.reshape(-1, self.boxvec.size)
        ave2 = dx.sum(0) / (dx.shape[0])
        self.transform.translate(x2, ave2)

        dist, x2 = self.finalize_best_match(x1, x2)    
        return dist, x1.ravel(), x2.ravel()  
    
    def __call__(self, coords1, coords2): 
        return self.align_fragments(coords1, coords2)    

    def finalize_best_match(self, x1, best_x2):
        ''' do final processing of the best match '''
        # Calculate the periodic distance between the two structures
        dist = self.measure.get_dist(x1, best_x2)
        return dist, best_x2.ravel()

