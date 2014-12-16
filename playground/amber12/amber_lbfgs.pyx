"""This file will define the interface for Amber12 with Cuda and the corresponding cuda lbfgs
"""

import numpy as np

from pele.potentials import BasePotential
from pele.potentials._pele cimport shared_ptr

#TODO: do these properly.  We need the iso c functions for amber12 cuda and the stuff for lbfgs
#cdef extern from "_mypotential.c":
cdef extern void mypot_energy_gradient(double * coords, int * natoms, double * e,
                                       double * grad, double *eps, double *sig )

# TODO: these need fixing.
# note: we should probably put all this stuff in a wales namespace
cdef extern from "filename":
    cdef cppclass c_gpu_amber "gpu_amber":
        gpu_amber(int ndof) 
cdef extern from "file.hpp" namespace "gpu_lbfgs":
    cdef cppclass  c_cuda_lbfgs "gpu_lbfgs::lbfgs":
        c_cuda_lbfgs(c_gpu_amber & potential, ...) except +
        minimize() except +


class Amber12CudaPotential(BasePotential):
    """this provides direct access to amber 12 and allows us to initialize it
    
    warning: probably multiple instantiations of this will refer to the same underlying Amber memory.
    Also, we need to ensure that the cuda LBFGS defined below accesses the memory prepared by this
    potential.
    
    This does not need to be a cdef class (I think) because we don't need to store any pure c objects
    """
    
    def __init__(self, natoms):
        self.natoms = natoms
        self._initialize_amber()
    
    def _initialize_amber(self):
        raise NotImplementedError()
    
    def getEnergyGradient(self, np.ndarray[double, ndim=1] coords):
        cdef np.ndarray[double, ndim=1] grad = np.zeros(coords.size)
#        cdef double * gdata = grad.data
#        cdef double * xdata = coords.data
        cdef double energy
        cdef double eps = 1.
        cdef double sig = 1.
        cdef int natoms = coords.size / 3
        mypot_energy_gradient(<double*> coords.data, <int *> &natoms, <double *> &energy, 
                    <double *>grad.data, <double*> &eps, <double*> &sig)
        return energy, grad
    
    def getEnergy(self, coords):
        e, g = self.getEnergyGradient(coords)
        return e
        

cdef class _cdef_LBFGSAmber12Cuda(object):
    """this provides access to the cuda lbfgs for amber
    
    This is specific to amber. perhaps it can be generalized later
    
    The constructor of the underlying cuda LBFGS accepts any object that
    derives from cost_function.  There is already such an object defined
    for amber12 called gpu_amber
    
    So first we have to create an object of type gpu_amber.  We then pass it
    to gpu_lbfgs.
    
    note: this needs to implement the same function as the other optimizers, e.g. run, get_result, etc.
    """
    cdef shared_ptr[gpu_amber] amber_cost_function # these are saved so the memory is not deleted
    cdef shared_ptr[gpu_lbfgs_lbfgs] amber_lbfgs
    cdef np.ndarray[double, ndim=1] coords
    cdef __cinit__(self, coords):
        self.coords = coords.copy()
        ndof = coords.size()
        # allocate memory for a new gpu_amber object and instantiate it.
        cdef gpu_amber * aptr = new gpu_amber(ndof)
        # the shared ptr will take ownership of the memory.  So you don't have to worry about deleting it.
        self.amber_cost_function = shared_ptr[gpu_amber](aptr) 
        
        # self.amber_cost_function.get() # returns a pointer to the actual amber cost function object
        cdef c_gpu_lbfgs * lbfgs_ptr = new c_gpu_lbfgs(* self.amber_cost_function.get())
        self.amber_lbfgs = shared_ptr[c_gpu_lbfgs](lbfgs_ptr)
    
    cdef run(self):
        self.amber_lbfgs.get().minimize(self.coords.data())
        # minimize will change the data in coords
        
        return self.get_result()
    
    def get_result():
        raise NotImplementedError()


class LBFGSAmber12Cuda(_cdef_LBFGSAmber12Cuda):
    # this is a python object and is easier to work with than a cdef class
