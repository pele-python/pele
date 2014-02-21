#cython: boundscheck=False
#cython: wraparound=False
##aaacython: noncheck=True

import numpy as np
cimport numpy as np 
from pele.potentials import _pele, _pythonpotential
from pele.potentials cimport _pele
from pele.potentials import _pythonpotential
from libcpp cimport bool as cbool

#===============================================================================
# shared pointer
#===============================================================================
cdef extern from "<memory>" namespace "std":
    cdef cppclass shared_ptr[T]:
        shared_ptr(T*)
        # Note: operator->, operator= are not supported

#===============================================================================
# pele::TakeStep
#===============================================================================

cdef extern from "pele/mc.h" namespace "pele":
    cdef cppclass cppTakeStep "pele::TakeStep"

cdef class _Cdef_TakeStep(object):
    """This class is the python interface for the c++ pele::TakeStep base class implementation
    """
    cdef cppTakeStep *thisptr
    
#===============================================================================
# pele::AcceptTest
#===============================================================================

cdef extern from "pele/mc.h" namespace "pele":
    cdef cppclass cppAcceptTest "pele::AcceptTest"
        
cdef class _Cdef_AcceptTest(object):
    """This class is the python interface for the c++ pele::AcceptTest base class implementation
    """
    cdef cppAcceptTest *thisptr
        
#===============================================================================
# pele::ConfTest
#===============================================================================

cdef extern from "pele/mc.h" namespace "pele":
    cdef cppclass cppConfTest "pele::ConfTest"

cdef class _Cdef_ConfTest(object):
    """This class is the python interface for the c++ pele::ConfTest base class implementation
    """
    cdef cppConfTest *thisptr

#===============================================================================
# pele::Action
#===============================================================================

cdef extern from "pele/mc.h" namespace "pele":
    cdef cppclass cppAction "pele::Action"
        
cdef class _Cdef_Action(object):
    """This class is the python interface for the c++ pele::Action base class implementation
    """
    cdef cppAction *thisptr

#===============================================================================
# pele::MC
#===============================================================================

cdef extern from "pele/mc.h" namespace "pele":
    cdef cppclass cppMC "pele::MC":
        cppMC(_pele.cBasePotential *, _pele.Array[double]&, double, double) except +
        void one_iteration() except +
        void run(size_t) except +
        void set_temperature(double) except +
        void set_stepsize(double) except +
        void add_action(shared_ptr[cppAction]) except +
        void add_accept_test( shared_ptr[cppAcceptTest]) except +
        void add_conf_test( shared_ptr[cppConfTest]) except +
        void set_takestep( shared_ptr[cppTakeStep]) except +
        void set_coordinates(_pele.Array[double], double) except +
        _pele.Array[double] get_coordinates() except +
        double get_energy() except +
        double get_accepted_fraction() except +

cdef class _Cdef_BaseMC(object):
    """This class is the python interface for the c++ pele::MC base class implementation
    """
    cdef cppMC* thisptr 