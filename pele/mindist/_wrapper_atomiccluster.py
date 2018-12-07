from __future__ import absolute_import
from .minpermdist_stochastic import MinPermDistCluster
from .exact_match import ExactMatchCluster
from ._minpermdist_policies import TransformAtomicCluster, MeasureAtomicCluster

class MinPermDistAtomicCluster(MinPermDistCluster):
    """ minpermdist for atomic cluster (3 carthesian coordinates per site)

    Parameters
    ----------

    permlist : optional
        list of allowed permutations. If nothing is given, all atoms will be
        considered as permutable. For no permutations give an empty list []
    can_invert : bool, optional
        also test for inversion

    See also
    --------

    MinPermDistCluster

    """
    def __init__(self, permlist=None, can_invert=True, **kwargs):
        transform=TransformAtomicCluster(can_invert=can_invert)
        measure = MeasureAtomicCluster(permlist=permlist)
        
        MinPermDistCluster.__init__(self, transform=transform, measure=measure, **kwargs)
        
class ExactMatchAtomicCluster(ExactMatchCluster):
    """ minpermdist for atomic cluster (3 carthesian coordinates per site)

    Parameters
    ----------

    permlist : optional
        list of allowed permutations. If nothing is given, all atoms will be
        considered as permutable. For no permutations give an empty list []
    can_invert : bool, optional
        also test for inversion

    See also
    --------

    ExactMatchCluster

    """
    def __init__(self, permlist=None, can_invert=True, **kwargs):
        transform = TransformAtomicCluster(can_invert=can_invert)
        measure = MeasureAtomicCluster(permlist=permlist)
        
        ExactMatchCluster.__init__(self, transform=transform, measure=measure, **kwargs)

