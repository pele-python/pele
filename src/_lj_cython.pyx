cimport _pygmin

cdef extern:
    void ljenergy_gradient_(double *x, int *natoms, double *e, double *grad,
        double *eps, double *sig, int *periodic, double *boxl)

cdef double lj_grad(double *x, double *grad, int n, void *userdata):
    cdef double *data = <double*>userdata
    cdef double sigma=data[0]
    cdef double eps=data[1]
    cdef double boxl=-1.0
    cdef int natoms = n / 3
    cdef double e
    cdef int periodic = 0
    ljenergy_gradient_(x, &natoms, &e, grad,
                       &eps, &sigma, &periodic, &boxl)
    return e

cdef double lj_energy(double *x, int n, void *userdata):
    pass

cdef class LJ_cython(_pygmin.Potential):
    cdef double _data[2]
    
    def __cinit__(self, sigma=1.0, epsilon=1.0):
        self._data[0] = sigma
        self._data[1] = epsilon
        self.thisptr = <_pygmin.cPotential*>new _pygmin.cPotentialFunction(
                                           &lj_energy,
                                           &lj_grad,
                                           <void*>self._data)
