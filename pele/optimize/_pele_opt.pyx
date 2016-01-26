# distutils: language = c++
import sys

import numpy as np

from pele.potentials import _pele
from pele.optimize import Result

cimport numpy as np
from libcpp cimport bool as cbool
cimport cython
import copy

from pele.potentials cimport _pele

@cython.boundscheck(False)
cdef pele_array_to_np_array(_pele.Array[double] v):
    """copy a pele Array into a new numpy array"""
    cdef np.ndarray[double, ndim=1] vnew = np.zeros(v.size(), dtype=float)
    cdef int i
    cdef int N = vnew.size
    for i in xrange(N):
        vnew[i] = v[i]
    
    return vnew


cdef class GradientOptimizer(object):
    """this class defines the python interface for c++ gradient optimizers 
    
    Notes
    -----
    for direct access to the underlying c++ optimizer use self.thisptr
    """
    res = None

    def one_iteration(self):
        self.thisptr.get().one_iteration()
        self.res = self.get_result()
        for event in self.events:
            event(coords=self.res.coords, energy=self.res.energy, rms=self.res.rms)
        return self.res
        
    def run(self, niter=None):
        if not self.events:
            # if we don't need to call python events then we can
            # go just let the c++ optimizer do it's thing
            if niter is None:
                self.thisptr.get().run()
            else:
                self.thisptr.get().run(niter)
            self.res = self.get_result()
        else:
            # we need to call python events after each iteration.
            if niter is None:
                niter = self.get_maxiter() - self.get_niter()
    
            for i in xrange(niter):
                if self.stop_criterion_satisfied():
                    break
                self.res = self.one_iteration()

        return copy.deepcopy(self.res)
            
    def reset(self, coords):
        cdef np.ndarray[double, ndim=1] ccoords = np.array(coords, dtype=float)
        self.thisptr.get().reset(_pele.Array[double](<double*> ccoords.data, ccoords.size))
    
    def stop_criterion_satisfied(self):
        return bool(self.thisptr.get().stop_criterion_satisfied())

    def get_maxiter(self):
        return self.thisptr.get().get_maxiter()

    def get_niter(self):
        return self.thisptr.get().get_niter()
    
    def get_result(self):
        """return a results object"""
        res = Result()
        
        cdef _pele.Array[double] xi = self.thisptr.get().get_x()
        cdef _pele.Array[double] gi = self.thisptr.get().get_g()
        x = pele_array_to_np_array(xi)
        g = pele_array_to_np_array(gi)
        
        res.energy = self.thisptr.get().get_f()
        res.coords = x
        res.grad = g
        
        res.rms = self.thisptr.get().get_rms()
        res.nsteps = self.thisptr.get().get_niter()
        res.nfev = self.thisptr.get().get_nfev()
        res.success = bool(self.thisptr.get().success())
        return res