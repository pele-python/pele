import numpy as np
cimport numpy as np
from pele.potentials import _pele, _pythonpotential
from pele.potentials cimport _pele
from pele.optimize import Result
cimport cython
from pele.potentials import _pythonpotential
import sys

# import the externally defined ljbfgs implementation
cdef extern from "_lbfgs.h" namespace "LBFGS_ns":
    cdef cppclass cppLBFGS "LBFGS_ns::LBFGS":
        cppLBFGS(_pele.cBasePotential *, _pele.Array[double] &, double, int) except +

        void run() except *
        void one_iteration() except *
        double get_f()
        _pele.Array[double] get_x()
        _pele.Array[double] get_g()
        
        void set_H0(double)
        void set_tol(double)
        void set_maxstep(double)
        void set_max_f_rise(double)
        void set_max_iter(int)
        void set_iprint(int)
        void set_verbosity(int)

        double get_rms()
        double get_H0()
        int get_nfev()
        int get_niter()
        int success()


cdef class LBFGS_CPP(object):
    """This class is the python interface for the c++ LBFGS implementation
    """
    cdef cppLBFGS *thisptr
    cdef _pele.BasePotential pot # this is stored so that the memory is not freed
    
    def __cinit__(self, x0, potential, double tol=1e-4, int M=4, double maxstep=0.1, 
                  double maxErise=1e-4, double H0=0.1, int iprint=-1, 
                  int nsteps=10000, int verbosity=0):
        if not issubclass(potential.__class__, _pele.BasePotential):
            print "LBFGS_CPP: potential is not subclass of BasePotential; wrapping it.", potential
            potential = _pythonpotential.CppPotentialWrapper(potential)
        cdef _pele.BasePotential pot = potential
        cdef np.ndarray[double, ndim=1] x0c = np.array(x0, dtype=float)
        self.thisptr = <cppLBFGS*>new cppLBFGS(pot.thisptr, 
               _pele.Array[double](<double*> x0c.data, x0c.size),
               tol, M)
        opt = self.thisptr
        opt.set_H0(H0)
        opt.set_maxstep(maxstep)
        opt.set_max_f_rise(maxErise)
        opt.set_max_iter(nsteps)
        opt.set_verbosity(verbosity)
        opt.set_iprint(iprint)
        self.pot = pot # so that the memory is not freed
        
    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr
            self.thisptr = NULL
        
    def run(self):
        self.thisptr.run()
        return self.get_result()

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
        res.H0 = self.thisptr.get_H0()
        res.success = self.thisptr.success()
        return res
    
    def one_iteration(self):
        self.thisptr.one_iteration()

