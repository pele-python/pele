import numpy as np
import itertools

__all__ = ["find_best_permutation"] 

have_minperm = False
have_hungarian = False
have_munkres = False

try:
    import minperm
    have_minperm = True    
except ImportError:
    pass
try:
    import hungarian
    have_hungarian = True
except ImportError:
    pass
try:    
    import munkres
    have_munkres = True
except ImportError:
    pass


_findBestPermutationList = None
if have_minperm:
    def _find_permutations(*args, **kwargs):
        return find_permutations_OPTIM(*args, **kwargs)

elif have_hungarian:
    def _find_permutations(*args, **kwargs):
        return find_permutations_hungarian(*args, **kwargs)
elif have_munkres:
    def _find_permutations(*args, **kwargs):
        return find_permutations_munkres(*args, **kwargs)
else:
    raise BaseException("No Hungarian algorithm implementation found!"
                        "Please compile minperm.f90 or install the hungarian or the munkres package")





def permuteArray(Xold, perm):
    #don't modify Xold
    Xnew = np.copy(Xold)
    permsorted = sorted(perm)
    for (iold, inew) in itertools.izip(permsorted, perm):
        #print iold, "->", inew
        Xnew[inew*3:inew*3+3] = Xold[iold*3:iold*3+3]

    return Xnew


def _make_cost_matrix(X1, X2):
    """
    return the cost matrix for use in the hungarian algorithm.
    
    the cost matrix is the distance matrix (squared) for all atoms in atomlist
    """
    cost = (((X1[np.newaxis,:] - X2[:,np.newaxis,:])**2).sum(2))
    return cost

def find_permutations_munkres( X1, X2, make_cost_matrix=_make_cost_matrix ):
    """
    For a given set of positions X1 and X2, find the best permutation of the
    atoms in X2.
    
    The positions must already be reshaped to reflect the dimensionality of the system!

    Use an implementation of the Hungarian Algorithm in the Python package
    index (PyPi) called munkres (another name for the algorithm).  The
    hungarian algorithm time scales as O(n^3), much faster than the O(n!) from
    looping through all permutations.

    http://en.wikipedia.org/wiki/Hungarian_algorithm
    http://pypi.python.org/pypi/munkres/1.0.5.2
    
    another package, hungarian, implements the same routine in comiled C
    http://pypi.python.org/pypi/hungarian/
    When I first downloaded this package I got segfaults.  The problem for me
    was casing an integer pointer as (npy_intp *).  I may add the corrected 
    version to pygmin at some point
    """
    #########################################
    # create the cost matrix
    # cost[j,i] = (X1(i,:) - X2(j,:))**2
    #########################################
    cost = make_cost_matrix(X1, X2)
    #cost = np.sqrt(cost)

    #########################################
    # run the munkres algorithm
    #########################################
    matrix = cost.tolist()
    m = munkres.Munkres()
    newind = m.compute(matrix)    
    
    #########################################
    # apply the permutation
    #########################################
    costnew = 0.
    X2new = np.copy(X2)
    new_indices = range(len(X1))
    for (iold, inew) in newind:
        costnew += cost[iold, inew]
        new_indices[inew] = iold

    dist = np.sqrt(costnew)
    return dist, new_indices

def find_permutations_hungarian( X1, X2, make_cost_matrix=_make_cost_matrix ):
    """
    For a given set of positions X1 and X2, find the best permutation of the
    atoms in X2.

    The positions must already be reshaped to reflect the dimensionality of the system!

    Use an implementation of the Hungarian Algorithm in the Python package
    index (PyPi) called munkres (another name for the algorithm).  The
    hungarian algorithm time scales as O(n^3), much faster than the O(n!) from
    looping through all permutations.

    http://en.wikipedia.org/wiki/Hungarian_algorithm
    http://pypi.python.org/pypi/munkres/1.0.5.2
    
    another package, hungarian, implements the same routine in comiled C
    http://pypi.python.org/pypi/hungarian/
    When I first downloaded this package I got segfaults.  The problem for me
    was casing an integer pointer as (npy_intp *).  I may add the corrected 
    version to pygmin at some point
    """
    #########################################
    # create the cost matrix
    # cost[j,i] = (X1(i,:) - X2(j,:))**2
    #########################################
    cost = make_cost_matrix(X1, X2)
    #cost = np.sqrt(cost)

    #########################################
    # run the hungarian algorithm
    #########################################
    newind1 = hungarian.lap(cost)
    perm = newind1[1]

    #note: the hungarian algorithm changes
    #the cost matrix.  I'm not sure why, and it may be a bug, 
    #but the indices it returns are still correct
#    if not np.all(cost >= 0):
#        m = np.max(np.abs(cost-costsave))
#        print "after hungarian cost greater than zero:, %g" % m
    
    
    #########################################
    # apply the permutation
    #########################################
    # TODO: how to get new distance?
    dist = -1
    return dist, perm

def find_permutations_OPTIM(X1, X2, boxl=None, make_cost_matrix=None):
    """
    use OPTIM's minperm() routine to calculate the optimum permutation
    """    
    
    if make_cost_matrix is not None:
        raise RuntimeError("cannot use a custom cost matrix with findBestPermutationListOPTIM")

    #deal with periodic boundary conditions
    periodic = boxl is not None
    if not periodic:
        #it must have a value for passing to fortran 
        boxl = 1.
    sx = sy = sz = boxl
        
    #run the minperm algorithm
    perm, dist, worstdist, worstradius = minperm.minperm(X1.flatten(), X2.flatten(), sx, sy, sz, periodic)
    perm -= 1 #fortran indexing

    #note, dist returned by minperm comes will only be accurate to 6 decimal places at best.
    #if we want a more accurate distance we should calculate it from the coordinates

    dist = np.sqrt(dist)
    return dist, perm


def find_best_permutation( X1, X2, permlist = None, user_algorithm=None, reshape=True, user_cost_matrix=None):
    """
    find the permutation of the atoms which minimizes the distance |X1-X2|
    
    With all the default parameters, findBestPermutation assumes that X1, X2
    are arrays of atoms in 3d space and performs reshaping on the coordinates. However,
    if you want to pass a 2D system or a custom array with own cost function, you can turn
    automatic reshaping off. 
    
    Parameters
    ----------
    X1, X2 : 
        the structures to align
    permlist : a list of lists
        A list of lists of atoms which are interchangable.
        e.g. for a 50/50 binary mixture, 
        
            permlist = [range(1,natoms/2), range(natoms/2,natoms)]
    user_algoriithm : None or callable
        you can optionally pass which algorithm to use.
    
    gen_cost_matrix : None or callable
        user function to generate the cost matrix
        
    reshape : boolean
        shall coordinate reshaping be performed.
    
    Returns
    -------
    dist : float
        the minimum distance
    perm:
        a list of all permutations
    
    Notes
    -----
    This for each list of interchangeable atoms in permlist the permutation
    which minimizes the distance between the two structures is found.  This minimimization
    is done by mapping the problem onto the linear assignment problem which can then be solved
    using graph theoretic techniques.  
    
    http://en.wikipedia.org/wiki/Linear_assignment_problem
    http://en.wikipedia.org/wiki/Hungarian_algorithm

    there are several packages in pypi which solve the linear assignment problem
    
    hungarian : c++ code wrapped in python.  scales roughly like natoms**2.5
    
    munkres : completely in python. scales roughly like natoms**3.  very slow for natoms > 10
    
    in addition we have wrapped the OPTIM version for use in pygmin.  It uses the sparse 
    version of the Jonker-Volgenant algorithm.  Furthermore the cost matrix calculated in 
    a compiled language for an additional speed boost. It scales roughtly like natoms**2

    """
    if reshape:
        X1 = X1.reshape([-1,3])
        X2 = X2.reshape([-1,3])
    
    if permlist is None:
        permlist = [range(len(X1))]
    
    newperm = range(len(X1))
    
    for atomlist in permlist:
        if user_algorithm is None:
            dist, perm = _find_permutations( X1[atomlist], X2[atomlist], make_cost_matrix=user_cost_matrix)
        else:
            dist, perm = user_algorithm( X1[atomlist], X2[atomlist], make_cost_matrix=user_cost_matrix)
            
        for atom,i in zip(atomlist,xrange(len(atomlist))):
            newperm[atom] = perm[i]
    return dist, newperm

