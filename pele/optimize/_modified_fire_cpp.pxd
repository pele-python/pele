from libcpp cimport bool as cbool
cimport numpy as np
cimport pele.optimize._pele_opt as _pele_opt
from pele.potentials cimport _pele
from pele.potentials._pele cimport shared_ptr

# import the externally defined modified_fire implementation
cdef extern from "pele/modified_fire.h" namespace "pele":
    cdef cppclass cppMODIFIED_FIRE "pele::MODIFIED_FIRE":
        cppMODIFIED_FIRE(shared_ptr[_pele.cBasePotential] , _pele.Array[double], 
                         double, double, double, size_t , double, double, 
                         double, double, double, cbool) except +