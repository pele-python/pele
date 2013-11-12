cimport pele.potentials._pele as _pele
import numpy as np
cimport numpy as np
from cpython cimport bool
from ctypes import c_size_t as size_t

# use external c++ class
cdef extern from "morse.h" namespace "pele":
    cdef cppclass  cMorse "pele::Morse":
        cMorse(double rho, double r0, double A) except +


cdef class Morse(_pele.BasePotential):
    """define the python interface to the c++ Morse implementation
    """
    def __cinit__(self, rho=1.0, r0=1.0, A=1.0):
        self.thisptr = <_pele.cBasePotential*>new cMorse(rho, r0, A)
