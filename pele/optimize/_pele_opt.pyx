import numpy as np
cimport numpy as np
from pele.potentials import _pele
from pele.potentials cimport _pele
from libcpp cimport bool as cbool
cimport cython
from pele.optimize import Result
import sys
# basic gradient optimizer interface stuff

# This is the base class for all gradient optimizers
cdef class GradientOptimizer(object):
    """this class defines the python interface for c++ gradient optimizers 
    
    Notes
    -----
    for direct access to the underlying c++ optimizer use self.thisptr
    """
    def __cinit__(self, *args, **kwargs):
        # store an instance to the current c++ class, will be used in every call
        self.thisptr = NULL #<cBasePotential*>new cBasePotential()
    
    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr
            self.thisptr = NULL
    
    def one_iteration(self):
        self.thisptr.one_iteration()
        
    def run(self, niter=None):
        if niter is None:
            self.thisptr.run()
        else:
            self.thisptr.run(niter)
        return self.get_result()
            
    def reset(self, coords):
        cdef np.ndarray[double, ndim=1] ccoords = np.array(coords, dtype=float)
        self.thisptr.reset(_pele.Array[double](<double*> ccoords.data, ccoords.size))
    
    def stop_criterion_satisfied(self):
        return bool(self.thisptr.stop_criterion_satisfied())

    def get_maxiter(self):
        return self.thisptr.get_maxiter()

    def get_niter(self):
        return self.thisptr.get_niter()
    
    @cython.boundscheck(False)
    def get_result(self):
        """return a results object"""
        res = Result()
        
        cdef _pele.Array[double] xi = self.thisptr.get_x()
        cdef _pele.Array[double] gi = self.thisptr.get_g()
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
        
        res.energy = self.thisptr.get_f()
        res.coords = x
        res.grad = g
        
        res.rms = self.thisptr.get_rms()        
        res.nsteps = self.thisptr.get_niter()
        res.nfev = self.thisptr.get_nfev()
        res.success = bool(self.thisptr.success())
        return res
