import numpy as np
cimport numpy as np
from playground.native_code import _pele
cimport _pele
from pele.optimize import Result

cdef extern from "_lbfgs.h" namespace "LBFGS_ns":
    cdef cppclass cppLBFGS "LBFGS_ns::LBFGS":
        cppLBFGS(_pele.cPotential *, _pele.Array &, double, int) except +

        void run() except +
        double get_f()
        _pele.Array get_x()
        _pele.Array get_g()
        
        void set_H0(double)
        void set_tol(double)
        void set_maxstep(double)
        void set_max_f_rise(double)
        void set_max_iter(int)
        void set_iprint(int)
        void set_verbosity(int)

        double get_f()
        double get_rms()
        double get_H0()
        int get_nfev()
        int get_niter()
        int success()


# we just need to set a different c++ class instance
cdef class LBFGS(object):
    cdef cppLBFGS *thisptr
    cdef _pele.Potential pot
    
    def __cinit__(self, _pele.Potential pot, np.ndarray[double, ndim=1, mode="c"] x0, tol=1e-4, int M=4,
                   double maxstep=0.1, double maxErise=1e-4, 
                   double H0=0.1, int iprint=-1, int nsteps=10000, int verbosity=0):
        self.thisptr = <cppLBFGS*>new cppLBFGS(pot.thisptr, 
                                               _pele.Array(<double*> x0.data, x0.size),
                                               tol,
                                               M)
        opt = self.thisptr
        opt.set_H0(H0)
        opt.set_maxstep(maxstep)
        opt.set_max_f_rise(maxErise)
        opt.set_max_iter(nsteps)
        opt.set_verbosity(verbosity)
        opt.set_iprint(iprint)
        self.pot = pot
        
    def __dealloc__(self):
        del self.thisptr
        
    def run(self):
        self.thisptr.run()
                       
#       x = np.ctypeslib.as_array(xptr, size=N)
#        if ret != 0 and ret != -1001: # ignore rounding errors
#            raise RuntimeError("lbfgs failed with error code %d:"%self.thisptr.error_code(), self.thisptr.error_string())    
        return self.get_result()

    def get_result(self):
        """return a results object"""
        res = Result()
        
        cdef _pele.Array xi = self.thisptr.get_x()
        cdef double *data = xi.data()
        cdef np.ndarray[double, ndim=1, mode="c"] x = np.zeros(xi.size())
        for i in xrange(xi.size()):
            x[i] = data[i]
        
        e, g = self.pot.get_energy_gradient(x)

        res.energy = e
        res.coords = x
        res.grad = g
        
        res.rms = self.thisptr.get_rms()        
        res.nsteps = self.thisptr.get_niter()
        res.nfev = self.thisptr.get_nfev()
        res.H0 = self.thisptr.get_H0()
        res.success = self.thisptr.success()
        return res

