from __future__ import absolute_import
from .aamindist import MeasureAngleAxisCluster
from pele.mindist.periodic_exact_match import TransformPeriodic
from pele.utils.rbtools import CoordsAdapter    
    
class MeasurePeriodicRigid(MeasureAngleAxisCluster): 
    """ Measure rules for rigid body systems with periodic boundary conditions

    Notes
    -----
    Ensure that the topology is of type RBTopologyBulk
    There is no meaningful definition of the centre of mass in a periodic
    system, so this is removed.
    Rotations are not yet implemented, so technically we should only use 
    non-cubic boxes.

    Parameters
    ----------
    topology: RBTopologyBulk
        the topology object for the system, containing the box size vector
    transform: Transform object, optional
        should be TransformPeriodicRigid().    
    permlist : list of lists, optional
        list of lists of identical atoms
    """      
    def get_com(self, X):
        raise NotImplementedError
    
    def find_rotation(self, X1, X2):
        raise NotImplementedError
    
    
class TransformPeriodicRigid(TransformPeriodic):
    ''' Interface for possible transformations on a set of rigid body 
    coordinates with periodic boundary conditions.
    
    The transform policy tells minpermdist how to perform translations only,
    because rotations and inversions are not generally relevant in periodic
    systems. Permutations should be implemented but currently are not.
    
    All transformations act in place, that means they change the current
    coordinates and do not make a copy.
    
    '''
    def translate(self,X,d):
        """ Performs an in-place translation of coordinates X by vector d """
        ca = CoordsAdapter(coords=X)
        if(ca.nrigid > 0):
            ca.posRigid += d
        if(ca.natoms > 0):
            ca.posAtom += d


class ExactMatchRigidPeriodic(object):
    """ Tests whether two rigid body periodic systems are identical
    
    Performs a translation to align the first atom in each set of
    coordinates, then tests whether the two structures are identical.
    
    Parameters
    ----------
    measure: Measure object
        should be MeasurePeriodicRigid(system.topology,TransformPeriodicRigid()
    accuracy: float (optional)
        the tolerance to which the structure comparison will test
    """      
    def __init__(self, measure, accuracy=.01):
        self.transform = TransformPeriodicRigid()
        self.measure = measure
        self.accuracy = accuracy        

    def __call__(self, x1, x2):      
        """ Tests whether coordinate sets x1 and x2 are translationally equivalent"""
        translate = x1[0:3]-x2[0:3]
        self.transform.translate(x2, translate)
#         print x1, x2
        dist = self.measure.get_dist(x1, x2)
        if dist <= self.accuracy:
            return True
        else:
            return False

