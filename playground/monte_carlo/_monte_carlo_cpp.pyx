import numpy as np
cimport numpy as np
from pele.potentials import _pele, _pythonpotential
from pele.potentials cimport _pele
from pele.optimize import Result
cimport cython
from pele.potentials import _pythonpotential
import sys
from libcpp cimport bool as cbool
from _pele_mc cimport *

cdef class _Cdef_MC(_Cdef_BaseMC):

    def __cinit__(self, _pele.BasePotential potential, coords, temperature, stepsize):
        cdef np.ndarray[double, ndim=1] coordsc = np.array(coords, dtype=float)        
        self.thisptr = <cppMC*>new cppMC(potential.thisptr, _pele.Array[double](<double*> coordsc.data, coordsc.size),
                                                                   temperature, stepsize)
    
    def add_action(self, _Cdef_Action action):
        self.thisptr.add_action(shared_ptr[cppAction](action.thisptr))
    
    def add_accept_test(self, _Cdef_AcceptTest test):
        self.thisptr.add_accept_test(shared_ptr[cppAcceptTest](test.thisptr))
    
    def add_conf_test(self, _Cdef_ConfTest test):
        self.thisptr.add_conf_test(shared_ptr[cppConfTest](test.thisptr))
    
    def set_takestep(self, _Cdef_TakeStep takestep):
        self.thisptr.set_takestep(shared_ptr[cppTakeStep](takestep.thisptr))
        
    def set_coordinates(self, coords, energy):
        cdef np.ndarray[double, ndim=1] ccoords = np.array(coords, dtype=float)
        self.thisptr.set_coordinates(_pele.Array[double](<double*> ccoords.data, ccoords.size), energy)
    
    def set_temperature(self, T):
        self.thisptr.set_temperature(T)
    
    def get_energy(self):
        energy = self.thisptr.get_energy()
        return energy
    
    @cython.boundscheck(False)
    def get_coords(self):
        """return a histogram array"""
        cdef _pele.Array[double] xi = self.thisptr.get_coords()
        cdef double *xdata = xi.data()
        cdef np.ndarray[double, ndim=1, mode="c"] x = np.zeros(xi.size())
        cdef size_t i
        for i in xrange(xi.size()):
            x[i] = xdata[i]
              
        return x
    
    def get_accepted_fraction(self):
        accepted_frac = self.thisptr.get_accepted_fraction()
        return accepted_frac
    
    def get_iterations_count(self):
        n = self.thisptr.get_iterations_count()
        return n
    
    def get_stepsize(self):
        s = self.thisptr.get_stepsize()
        return s
    
    def one_iteration(self):
        self.thisptr.one_iteration()
    
    def run(self, niter):
        self.thisptr.run(niter)

class MC(_Cdef_MC):
    """This class is the python interface for the c++ MC implementation.
    """