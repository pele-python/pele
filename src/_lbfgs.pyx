#import numpy as np
cimport numpy as np
import _pygmin
cimport _pygmin


cdef extern from "lbfgs_wrapper.cpp" namespace "pygmin":
    cdef cppclass cLBFGS "pygmin::LBFGS":
        cLBFGS() except +
        cLBFGS(_pygmin.cPotential *) except +
        #
        int run(_pygmin.Array &x)

# we just need to set a different c++ class instance
cdef class LBFGS:
    cdef cLBFGS *thisptr
    cdef _pygmin.Potential pot
    
    def __cinit__(self, _pygmin.Potential pot):
        self.thisptr = <cLBFGS*>new cLBFGS(pot.thisptr)
        self.pot = pot
        
    def __dealloc__(self):
        del self.thisptr
        
    def run(self, np.ndarray[double, ndim=1, mode="c"] x not None):
        ret = self.thisptr.run(_pygmin.Array(<double*> x.data, x.size))
        return x, self.pot.get_energy(x)
        