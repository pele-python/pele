"""
# distutils: language = C++
"""
import numpy as np

cimport numpy as np

cimport pele.potentials._pele as _pele
from pele.potentials._pele cimport shared_ptr

# use external c++ class
cdef extern from "pele/morse.h" namespace "pele":
    cdef cppclass  cMorse "pele::Morse":
        cMorse(double rho, double r0, double A) except +


cdef class Morse(_pele.SimplePairwisePotential):
    """define the python interface to the c++ Morse implementation
    """
    def __cinit__(self, rho=1.0, r0=1.0, A=1.0):
        self.thisptr = shared_ptr[_pele.cBasePotential]( <_pele.cBasePotential*>new cMorse(rho, r0, A) )
        self.spp_ptr = shared_ptr[_pele.cppSimplePairwisePotentialInterface](<_pele.cppSimplePairwisePotentialInterface*> self.thisptr.get())
