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
