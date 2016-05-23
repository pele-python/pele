from libcpp cimport bool as cbool

cimport pele.optimize._pele_opt as _pele_opt
from pele.potentials cimport _pele
from pele.potentials._pele cimport shared_ptr

cdef extern from "pele/stochastic_gradient_descent.h" namespace "pele":
    cdef cppclass cppStochasticGradientDescent "pele::StochasticGradientDescent":
        cppStochasticGradientDescent(shared_ptr[_pele.cBasePotential],
            _pele.Array[double], double, double, size_t, cbool) except+
