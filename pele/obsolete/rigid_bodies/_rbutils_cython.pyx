import numpy as np
cimport numpy as np
cimport cython


@cython.boundscheck(False)
def _moleculeGetGradients(np.ndarray[np.float64_t, ndim=3] drmat,
                           np.ndarray[np.float64_t, ndim=2] sitexyz_molframe, 
                           np.ndarray[np.float64_t, ndim=2] sitegrad,
                           Py_ssize_t nsites):
    cdef np.ndarray[np.float64_t, ndim=1] aagrad = np.zeros(3, dtype=np.float64)
    cdef Py_ssize_t i, j, k, isite
    for k in range(3):
        for i in range(3):
            for j in range(3):
                for isite in range(nsites):
                    aagrad[k] += drmat[k,i,j] * sitexyz_molframe[isite,j] * sitegrad[isite,i]

    return aagrad
