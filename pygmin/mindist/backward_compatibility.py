from minpermdist_stochastic import MinPermDistCluster
from _minpermdist_policies import MeasureAtomicCluster
from permutational_alignment import optimize_permutations

def minPermDistStochastic(X1, X2, niter=100, permlist=None, verbose=False, accuracy=0.01,
                      check_inversion=True):
    if check_inversion is False:
        raise NotImplementedError
    
    print "WARNING: minPermDistStochastic is obsolete"
    
    return MinPermDistCluster(niter=niter, measure=MeasureAtomicCluster(permlist=permlist),
                             accuracy = accuracy, verbose=verbose)(X1, X2)
    
def findBestPermutation( X1, X2, permlist = None, user_algorithm=None):
    print "WARNING: findBestPermutation is obsolete, use optimize_permutations instead"    
    return optimize_permutations(X1, X2, permlist=permlist, user_algorithm=user_algorithm)
    return dist, X1, X2new.flatten()

class MinDistWrapper(object):
    """
    wrap a mindist routine into a callable object with the form mindist(X1, X2)
    
    Parameters
    ----------
    mindist : callable
        the mindist routine
    args : 
        extra arguements for mindist
    kwargs : 
        extra keyword arguments for mindist
    """
    def __init__(self, mindist, *args, **kwargs):
        self.mindist = mindist
        self.args = args
        self.kwargs = kwargs
    
    def __call__(self, X1, X2):
        return self.mindist(X1, X2, *self.args, **self.kwargs)