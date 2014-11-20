"""
This is the preferred method for implementing a potential in pele.  It connects
to the c++ backend that parts of pele is built on.  This means that the
implementation is very straightforward, but the most important benefit of doing
it this way is that it can be much faster.  With this method the c++ potential
is passed directly to the c++ LBFGS routine and the entire optimization is done
in c++ with no need for any slow python calls.  

Compilation
-----------
This cython file must be cythonized.  We must pass as an include directory the
place where the cython wrapped c++ BasePotential is defined
    
    cython --cplus mypotential.pyx -I ../../../pele/potentials/

Now we compile mypotential.cpp and _mypotential.hpp into a shared object
library which can be called from python
"""
cimport cython
#import pele.potentials._pele as _pele
cimport pele.potentials._pele as _pele
from pele.potentials._pele cimport BasePotential
# note: it is required to explicitly import BasePotential.  The compilation
# will fail if you try to use it as _pele.BasePotential.  I don't know why this
# is

# use external c++ class
cdef extern from "_mypotential.hpp" namespace "pele":
    cdef cppclass  cMyPot "pele::MyPot":
        cMyPot(double sig, double eps) except +

cdef class MyPotCpp(BasePotential):
    """define the python interface to the c++ MyPot implementation
    """
    def __cinit__(self, natoms, sig=1.0, eps=1.0):
        self.thisptr = _pele.shared_ptr[_pele.cBasePotential](<_pele.cBasePotential*> new cMyPot(sig, eps))
