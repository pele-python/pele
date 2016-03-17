cimport pele.optimize._pele_opt as _pele_opt
from pele.potentials cimport _pele
from pele.potentials._pele cimport shared_ptr

cdef extern from "pele/steepest_descent.h" namespace "pele":
    cdef cppclass cppSteepestDescent "pele::SteepestDescent":
        cppSteepestDescent(shared_ptr[_pele.cBasePotential],
            _pele.Array[double], double, double) except+
