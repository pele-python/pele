from playground.native_code cimport _pygmin

# use external c++ class
cdef extern from "lj.h" namespace "pygmin":
    cdef cppclass  cLJ "pygmin::LJ":
        cLJ(double C6, double C12) except +

# we just need to set a different c++ class instance
cdef class LJ(_pygmin.Potential):
    def __cinit__(self, eps=1.0, sigma=1.0):
        self.thisptr = <_pygmin.cPotential*>new cLJ(4.*eps*sigma**6, 4.*eps*sigma**12)
