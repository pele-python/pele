from libcpp cimport bool as cbool

cimport pele.optimize._pele_opt as _pele_opt
from pele.potentials cimport _pele
from pele.potentials._pele cimport shared_ptr

cdef extern from "pele/stochastic_diagonal_levenberg_marquardt.h" namespace "pele":
    cdef cppclass cppStochasticDiagonalLevenbergMarquardt "pele::StochasticDiagonalLevenbergMarquardt":
        cppStochasticDiagonalLevenbergMarquardt(
            shared_ptr[_pele.cBasePotential], _pele.Array[double],
            double, double, double, double, size_t, cbool) except +
