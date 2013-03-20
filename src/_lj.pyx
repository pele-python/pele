cimport _pygmin

# use external c++ class
cdef extern from "lj.h" namespace "pygmin":
    cdef cppclass  cLJ "pygmin::LJ":
        cLJ(double C6, double C12) except +

# we just need to set a different c++ class instance
cdef class LJ(_pygmin.Potential):
    def __cinit__(self, C6=1.0, C12=1.0):
        self.thisptr = <_pygmin.cPotential*>new cLJ(C6, C12)
