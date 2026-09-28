# cython: language_level=3str
# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
from libcpp cimport bool as cbool
cimport pele.optimize._pele_opt as _pele_opt
from pele.potentials cimport _pele
from pele.potentials._pele cimport shared_ptr
from libcpp.vector cimport vector

cdef extern from "pele/optimizer.hpp" namespace "pele":
    cdef enum StopCriterionType:
        GRADIENT,
        STEPNORM,
        NEWTON

# import the externally defined modified_fire implementation
cdef extern from "pele/modified_fire.hpp" namespace "pele":
    cdef cppclass cppMODIFIED_FIRE "pele::MODIFIED_FIRE":
        cppMODIFIED_FIRE(shared_ptr[_pele.cBasePotential] , _pele.Array[double], 
                         double, double, double, size_t , double, double, 
                         double, double, double, cbool, cbool, int, StopCriterionType stop_criterion_type) except +
        vector[double] get_time_trajectory() except +
        vector[double] get_gradient_norm_trajectory() except +
        vector[double] get_distance_trajectory() except +
        vector[double] get_energy_trajectory() except +
        vector[double] get_costly_time_trajectory() except +
        vector[vector[double]] get_coordinate_trajectory() except +
        vector[vector[double]] get_gradient_trajectory() except +