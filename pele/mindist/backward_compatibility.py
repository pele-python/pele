from __future__ import absolute_import
from .minpermdist_stochastic import MinPermDistCluster
from ._minpermdist_policies import MeasureAtomicCluster
from .permutational_alignment import optimize_permutations
from pele.utils import rotations
from pele.utils import rotations as rot
from .rmsfit import findrotation
import numpy as np

# def getAlignRotation(XA, XB):
#     print "WARNING: getAlignRotation is obsolete, use find_rotation"
#     return rotations.mx2q(findrotation(XA, XB))

# def minPermDistStochastic(X1, X2, niter=100, permlist=None, verbose=False, accuracy=0.01,
#                       check_inversion=True):
#     if check_inversion is False:
#         raise NotImplementedError
#     
#     print "WARNING: minPermDistStochastic is obsolete"
#     
#     return MinPermDistCluster(niter=niter, measure=MeasureAtomicCluster(permlist=permlist),
#                              accuracy = accuracy, verbose=verbose)(X1, X2)
    
# def findBestPermutation( X1, X2, permlist = None, user_algorithm=None):
#     print "WARNING: findBestPermutation is obsolete, use optimize_permutations instead"    
#     return optimize_permutations(X1, X2, permlist=permlist, user_algorithm=user_algorithm)
#     return dist, X1, X2new.flatten()

# def alignCoM( X1, X2):
#     """
#     align the center of mass of X2 with that of X1
#     """
#     natoms = len(X1) / 3
#     for i in xrange(3):
#         com = np.sum( X1[i::3] - X2[i::3] )
#         com /= natoms
#         X2[i::3] += com

def CoMToOrigin( X1):
    """
    move the center of mass to the origin
    """
    X1 = np.reshape(X1, [-1,3])
    natoms = X1.shape[0]
    com = X1.sum(0) / natoms
    X1 -= com
    return X1.reshape(-1)

# def alignRotation(XA, XB):
#     """
#     Align structure XB with XA
# 
#     Align structure XB to be as similar as possible to structure XA.
#     To be precise, rotate XB, so as to minimize the distance |XA - XB|.
# 
#     Rotations will be done around the origin, not the center of mass
#     """
#     nsites = len(XA)/3
# 
#     dist, Q2 = getAlignRotation(XA, XB)
#     ###################################################################
#     # Q2 contains the quaternion which rotates XB to best align with X1.
#     # rotate XB according to Q2
#     ###################################################################
# 
#     rot_mx = rot.q2mx( Q2 )
#     #print rot_mx
#     for j in range(nsites):
#         i = 3*j
#         XB[i:i+3] = np.dot( rot_mx, XB[i:i+3] )
#     
#     dist = np.linalg.norm(XA - XB) #more precise than what it returned
# 
#     return dist, XB

# def getDistxyz( xyz1, xyz2 ):
#     return np.linalg.norm(xyz1 - xyz2)

# class MinDistWrapper(object):
#     """
#     wrap a mindist routine into a callable object with the form mindist(X1, X2)
#     
#     Parameters
#     ----------
#     mindist : callable
#         the mindist routine
#     args : 
#         extra arguements for mindist
#     kwargs : 
#         extra keyword arguments for mindist
#     """
#     def __init__(self, mindist, *args, **kwargs):
#         self.mindist = mindist
#         self.args = args
#         self.kwargs = kwargs
#     
#     def __call__(self, X1, X2):
#         return self.mindist(X1, X2, *self.args, **self.kwargs)

