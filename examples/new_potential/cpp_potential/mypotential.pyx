"""
This is the preferred method for implementing a potential in pele.  It connects
to the c++ backend that parts of pele is built on.  This means that the
implementation is very straightforward, but the most important benefit of doing
it this way is that it can be much faster.  With this method the c++ potential
is passed directly to the c++ LBFGS routine and the entire optimization is done
in c++ with no need for any slow python calls.  
"""

cimport cython
cimport pele.potentials._pele as _pele
from pele.potentials._pele cimport BasePotential

cdef extern from "_mypotential.hpp" namespace "pele":
    cdef cppclass  cMyPot "pele::MyPot":
        cMyPot(double sig, double eps) except +

cdef class MyPotCpp(BasePotential):
    """define the python interface to the c++ MyPot implementation
    """
    def __cinit__(self, natoms, sig=1.0, eps=1.0):
        self.thisptr = _pele.shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new cMyPot(sig, eps))

