# This file demonstrates a pure cython wrapper for a fortran function
#
# I didn't figure out yet how to declare inherited c++ classes in python,
# therefore this approach uses the PotentialFunction class to add a
# potential callback. 
#
# Sigma and epsilon have to be passed as userdata since callback functions
# cannot have any instance data and any variable data has to be passed in
# the function call.

cimport _pygmin

# this is implemented in external fortran code
cdef extern:
    void ljenergy_gradient_(double *x, int *natoms, double *e, double *grad,
        double *eps, double *sig, int *periodic, double *boxl)

# we need to do some wrapping or arguments
# the wrapping can be simplified a lot with little modification to 
# fortran code (change to c type declaration)
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

# energy callback not yet implemented
cdef double lj_energy(double *x, int n, void *userdata):
    pass

# define the potential class
cdef class LJ_cython(_pygmin.Potential):
    cdef double _data[2]
    
    def __cinit__(self, sigma=1.0, epsilon=1.0):
        self._data[0] = sigma
        self._data[1] = epsilon
        # PotentialFunction uses lj_energy and lj_grad as callback
        # and self._data is passed to each function call
        self.thisptr = <_pygmin.cPotential*>new _pygmin.cPotentialFunction(
                                           &lj_energy,
                                           &lj_grad,
                                           <void*>self._data)
