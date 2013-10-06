import sys
import numpy as np
cimport numpy as np
cimport cython


@cython.boundscheck(False) # turn of bounds-checking for entire function
@cython.wraparound(False)
#@cython.
def _compute_LBFGS_step(
                           np.ndarray[double, ndim=1] G,
                           np.ndarray[double, ndim=2] s,
                           np.ndarray[double, ndim=2] y,
                           np.ndarray[double, ndim=1] rho,
                           long int k, double H0
                           ):
    cdef long int i, M, N, istart, j1, j2, nindices
    cdef double sq, yz, beta
    cdef np.ndarray[double, ndim=1] stp
    cdef np.ndarray[double, ndim=1] a
    cdef np.ndarray[long int, ndim=1] indices
    
    M = s.shape[0]
    N = G.size
    assert s.shape[0] == M
    assert s.shape[1] == N
#    assert s.shape == y.shape
    assert y.shape[1] == N
    assert rho.size == M
    
    stp = G.copy()
    a = np.zeros(M)

    if k == 0:
        #make first guess for the step length cautious
        gnorm = np.linalg.norm(G)
        stp *= -H0 * min(gnorm, 1. / gnorm)
        return stp
    
    indices = np.array([ i % M for i in range(max([0, k - M]), k, 1) ], np.integer)
    nindices = len(indices)
    
    
    # loop through the history, most recent first
    for j1 in xrange(nindices):
        i = indices[nindices - j1 - 1]
#        a[i] = rho[i] * np.dot( s[i,:], q )
#        q -= a[i] * y[i,:]
        sq = 0.
        for j2 in xrange(N):
            sq += s[i,j2] * stp[j2]
        a[i] = rho[i] * sq
        for j2 in xrange(N):
            stp[j2] -= a[i] * y[i,j2]
    
    # include our estimate for diagonal component of the inverse hessian
    for j2 in xrange(N):
        stp[j2] *= H0
    
    # loop through the history, most recent last
    for j1 in xrange(nindices):
        i = indices[j1]
#        beta = rho[i] * np.dot( y[i,:], z )
#        z += s[i,:] * (a[i] - beta)
        yz = 0.
        for j2 in xrange(N):
            yz += y[i,j2] * stp[j2]
        beta = rho[i] * yz
        for j2 in xrange(N):
            stp[j2] += s[i,j2] * (a[i] - beta)
    
        
    # step should point downhill
#    stp = - stp
    for j2 in xrange(N):
        stp[j2] = -stp[j2]


    return stp
    
