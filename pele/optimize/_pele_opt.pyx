# distutils: language = c++
import sys

import numpy as np

from pele.potentials import _pele
from pele.optimize import Result

cimport numpy as np
from libcpp cimport bool as cbool
cimport cython

from pele.potentials cimport _pele

cdef class GradientOptimizer(object):
    """this class defines the python interface for c++ gradient optimizers 
    
    Notes
    -----
    for direct access to the underlying c++ optimizer use self.thisptr
    """
    def one_iteration(self):
        self.thisptr.get().one_iteration()
        res = self.get_result()
        for event in self.events:
            event(coords=res.coords, energy=res.energy, rms=res.rms)
        return res
        
    def run(self, niter=None):
        if not self.events:
            # if we don't need to call python events then we can
            # go just let the c++ optimizer do it's thing
            if niter is None:
                self.thisptr.get().run()
            else:
                self.thisptr.get().run(niter)
        else:
            # we need to call python events after each iteration.
            if niter is None:
                niter = self.get_maxiter() - self.get_niter()
    
            for i in xrange(niter):
                if self.stop_criterion_satisfied():
                    break
                self.one_iteration()
            
        return self.get_result()
            
    def reset(self, coords):
        cdef np.ndarray[double, ndim=1] ccoords = np.array(coords, dtype=float)
        self.thisptr.get().reset(_pele.Array[double](<double*> ccoords.data, ccoords.size))
    
    def stop_criterion_satisfied(self):
        return bool(self.thisptr.get().stop_criterion_satisfied())

    def get_maxiter(self):
        return self.thisptr.get().get_maxiter()

    def get_niter(self):
        return self.thisptr.get().get_niter()
    
    @cython.boundscheck(False)
    def get_result(self):
        """return a results object"""
        res = Result()
        
        cdef _pele.Array[double] xi = self.thisptr.get().get_x()
        cdef _pele.Array[double] gi = self.thisptr.get().get_g()
        cdef double *xdata = xi.data()
        cdef double *gdata = gi.data()
        cdef np.ndarray[double, ndim=1, mode="c"] x = np.zeros(xi.size())
        cdef np.ndarray[double, ndim=1, mode="c"] g = np.zeros(xi.size())
        cdef size_t i
        for i in xrange(xi.size()):
            x[i] = xdata[i]
            g[i] = gdata[i]

        #jake> it's anoying having to copy the data manually like this.
        # We can possibly use np.frombuffer(), thought I haven't gotten it to work.
        # We can also maybe use the c function PyArray_SimpleNewFromData.
        # In the meantime the loop won't be too slow if we use cython properly to speed it up.
        #cdef np.ndarray[double, ndim=1, mode="c"] g2
        #g2 = np.frombuffer(xi.data(), dtype=np.float64, count=xi.size())
        
        res.energy = self.thisptr.get().get_f()
        res.coords = x
        res.grad = g
        
        res.rms = self.thisptr.get().get_rms()        
        res.nsteps = self.thisptr.get().get_niter()
        res.nfev = self.thisptr.get().get_nfev()
        res.success = bool(self.thisptr.get().success())
        return res
