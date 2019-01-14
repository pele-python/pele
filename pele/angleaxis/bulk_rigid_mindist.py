from __future__ import print_function
import numpy as np
from pele.mindist.periodic_exact_match import TransformPeriodic
from pele.utils.rbtools import CoordsAdapter
from inspect import stack

class MinDistBulkRigid(object):
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
                          
    def align_fragments(self, coords1, coords2): 
        """
        Obtain the best alignment between two configurations of a periodic system
        
        Parameters
        ----------
        coords1, coords2 : np.array 
            the structures to align.  X2 will be aligned with X1
            Both structures are arrays of rigid body com positions and aa vectors

        Returns
        -------
        a triple of (dist, coords1, coords2). coords1 are the unchanged coords1
        and coords2 are brought in best alignment with coords2
        """

        if self.verbose:
            print("Measure:")
            print(self.measure)
            print("Transform:")
            print(self.transform)
            print("Measure.topology:")
            print(self.measure.topology)
            print("Called by", stack())
        
        # we don't want to change the given coordinates
        coords1 = coords1.copy()
        coords2 = coords2.copy()

        x1 = np.copy(coords1)
        x2 = np.copy(coords2)

        self.distbest = self.measure.get_dist(x1, x2)
        ca1 = CoordsAdapter(coords=x1)
        ca2 = CoordsAdapter(coords=x2)              
          
        dx = ca1.posRigid - ca2.posRigid
        dx -= np.round(dx / self.boxvec) * self.boxvec
        ave2 = dx.sum(0) / ca1.nrigid 
        self.transform.translate(x2, ave2)

        dist, x2 = self.finalize_best_match(coords1, x2)    
        return dist, coords1, x2  
    
    def __call__(self, coords1, coords2): 
        return self.align_fragments(coords1, coords2)    

    def finalize_best_match(self, x1, best_x2):
        ''' do final processing of the best match '''
        ca1 = CoordsAdapter(coords=x1)     
        ca2 = CoordsAdapter(coords=best_x2)
        dx = ca1.posRigid - ca2.posRigid
        dx = np.round(dx / self.boxvec) * self.boxvec
        self.transform.translate(best_x2, dx)

        dist = self.measure.get_dist(x1, best_x2)
#         if (dist - self.distbest) > 1e-6:
#             raise RuntimeError(dist, self.distbest, "Permutational alignment has increased the distance metric")        
        if self.verbose:
            print("finaldist", dist, "distmin", self.distbest)
        return dist, best_x2

