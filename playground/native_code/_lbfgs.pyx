import numpy as np
cimport numpy as np
from playground.native_code import _pele
cimport _pele


cdef extern from "lbfgs_wrapper.cpp" namespace "pele":
    cdef cppclass cLBFGS "pele::LBFGS":
        cLBFGS() except +
        cLBFGS(_pele.cPotential *) except +
        #
        int run(_pele.Array &x)

        const char * error_string()
        int error_code()

# we just need to set a different c++ class instance
cdef class LBFGS:
    cdef cLBFGS *thisptr
    cdef _pele.Potential pot
    
    def __cinit__(self, _pele.Potential pot):
        self.thisptr = <cLBFGS*>new cLBFGS(pot.thisptr)
        self.pot = pot
        
    def __dealloc__(self):
        del self.thisptr
        
    def run(self, np.ndarray[double, ndim=1, mode="c"] x not None):
        ret = self.thisptr.run(_pele.Array(<double*> x.data, x.size))
        if ret != 0 and ret != -1001: # ignore rounding errors
            raise RuntimeError("lbfgs failed with error code %d:"%self.thisptr.error_code(), self.thisptr.error_string())
        e, g = self.pot.get_energy_gradient(x)
        return x, e, np.linalg.norm(g) / np.sqrt(g.size), -1
        
