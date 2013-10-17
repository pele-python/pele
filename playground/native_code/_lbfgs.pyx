import numpy as np
cimport numpy as np
from playground.native_code import _pele
cimport _pele

cdef extern from "_lbfgs.h" namespace "LBFGS_ns":
    cdef cppclass cppLBFGS "LBFGS_ns::LBFGS":
        cppLBFGS(_pele.cPotential *, _pele.Array &, int) except +

        void run() except +
        double get_f()
        double* get_x()
        double* get_g()
        int get_N()


# we just need to set a different c++ class instance
cdef class LBFGS(object):
    cdef cppLBFGS *thisptr
    cdef _pele.Potential pot
    
    def __cinit__(self, _pele.Potential pot, np.ndarray[double, ndim=1, mode="c"] x0, M=4):
        cdef int Mint = M
        self.thisptr = <cppLBFGS*>new cppLBFGS(pot.thisptr, 
                                               _pele.Array(<double*> x0.data, x0.size), 
                                               Mint)
        self.pot = pot
        
    def __dealloc__(self):
        del self.thisptr
        
    def run(self):
        self.thisptr.run()
        N = self.thisptr.get_N()
        cdef double* xi = self.thisptr.get_x()
        N = int(self.thisptr.get_N())
        x = np.zeros(N)
        for i in xrange(N):
            x[i] = xi[i]
        
#       x = np.ctypeslib.as_array(xptr, size=N)
#        if ret != 0 and ret != -1001: # ignore rounding errors
#            raise RuntimeError("lbfgs failed with error code %d:"%self.thisptr.error_code(), self.thisptr.error_string())
        e, g = self.pot.get_energy_gradient(x)
        return x, e, np.linalg.norm(g) / np.sqrt(g.size), -1
        
