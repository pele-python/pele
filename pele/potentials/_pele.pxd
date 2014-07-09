
#===============================================================================
# shared pointer
#===============================================================================
cdef extern from "<memory>" namespace "std":
    cdef cppclass shared_ptr[T]:
        shared_ptr() except+
        shared_ptr(T*) except+
        T* get() except+
        # T& operator*() # doesn't do anything
        # Note: operator->, operator= are not supported

#===============================================================================
# pele::Array
#===============================================================================
cdef extern from "pele/array.h" namespace "pele":
    cdef cppclass Array[dtype] :
        Array() except +
        Array(dtype*, size_t n) except +
        size_t size() except +
        dtype *data() except +
        dtype & operator[](size_t) except +

#===============================================================================
# pele::BasePotential
#===============================================================================
cdef extern from "pele/base_potential.h" namespace "pele":
    cdef cppclass  cBasePotential "pele::BasePotential":
        cBasePotential() except +
        double get_energy(Array[double] &x) except +
        double get_energy_gradient(Array[double] &x, Array[double] &grad) except +
        double get_energy_gradient_hessian(Array[double] &x, Array[double] &g, Array[double] &hess) except +
        void get_hessian(Array[double] &x, Array[double] &hess) except +
        void numerical_gradient(Array[double] &x, Array[double] &grad, double eps) except +
        void numerical_hessian(Array[double] &x, Array[double] &hess, double eps) except +

#cdef extern from "potentialfunction.h" namespace "pele":
#    cdef cppclass  cPotentialFunction "pele::PotentialFunction":
#        cPotentialFunction(
#            double (*energy)(Array[double] x, void *userdata) except *,
#            double (*energy_gradient)(Array[double] x, Array[double] grad, void *userdata) except *,
#            void *userdata) except +

#===============================================================================
# cython BasePotential
#===============================================================================
cdef class BasePotential:
    cdef shared_ptr[cBasePotential] thisptr      # hold a C++ instance which we're wrapping

#===============================================================================
# pele::CombinedPotential
#===============================================================================
cdef extern from "pele/combine_potentials.h" namespace "pele":
    cdef cppclass  cCombinedPotential "pele::CombinedPotential":
        cCombinedPotential() except +
        double get_energy(Array[double] &x) except +
        double get_energy_gradient(Array[double] &x, Array[double] &grad) except +
        void add_potential(shared_ptr[cBasePotential] potential) except +
    
