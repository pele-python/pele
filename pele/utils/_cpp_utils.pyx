"""
# distutils: language = C++
"""
import numpy as np

cimport numpy as np

# cython has no support for integer template argument.  This is a hack to get around it
# https://groups.google.com/forum/#!topic/cython-users/xAZxdCFw6Xs
# Basically you fool cython into thinking INT2 is the type integer,
# but in the generated c++ code you use 2 instead.
# The cython code MyClass[INT2] will create c++ code MyClass<2>.
cdef extern from *:
    ctypedef int INT3 "3"    # a fake type
    ctypedef int INT4 "4"   


cdef extern from "pele/vecn.h" namespace "pele":
    cdef cppclass VecN "pele::VecN"[N]:
        VecN()
        double & operator[](size_t) except +
    cdef cppclass Matrix33 "pele::MatrixNM<3,3>":
        Matrix33()
        double & operator()(size_t, i, size_t j) except +
        double * data() except +

ctypedef VecN[INT3] Vec3 
ctypedef VecN[INT4] Vec4 

# use external c++ class
cdef extern from "pele/rotations.h" namespace "pele":
    cdef Vec3 c_rotate_aa "pele::rotate_aa" (Vec3 &, Vec3 &) except +
    cdef Vec3 rot_mat_to_aa(Matrix33 & mx) except +
    cdef Vec4 aa_to_quaternion(Vec3 & mx) except +

cdef Vec3 to_vec(p):
    cdef Vec3 v
    assert len(p) == 3
    v[0] = p[0]
    v[1] = p[1]
    v[2] = p[2]
    return v

cdef np.ndarray[double] vec_to_np(Vec3 & p):
    v = np.zeros(3)
    v[0] = p[0]
    v[1] = p[1]
    v[2] = p[2]
    return v

cdef Vec4 to_vec4(p):
    cdef Vec4 v
    assert len(p) == 4
    v[0] = p[0]
    v[1] = p[1]
    v[2] = p[2]
    v[3] = p[3]
    return v

cdef np.ndarray[double] vec4_to_np(Vec4 & p):
    v = np.zeros(4)
    v[0] = p[0]
    v[1] = p[1]
    v[2] = p[2]
    v[3] = p[3]
    return v

cdef Matrix33 to_mat(A) except *:
    cdef Matrix33 Anew
    try:
        A = A.reshape(-1)
    except AttributeError:
        pass
    if len(A) != 9:
        print A
        raise ValueError("A should have length 9 but it is length " + str(len(A)) )
    cdef int i
    cdef double * data = Anew.data()
    for i in xrange(9):
        data[i] = A[i]
    return Anew

cdef np.ndarray[double] mat_to_np(Matrix33 & A):
    cdef np.ndarray[double, ndim=1] Anew = np.zeros(9)
    cdef int i
    cdef double * data = A.data()
    for i in xrange(9):
        Anew[i] = data[i]
    return Anew


cpdef rotate_aa(p1, p2):
    """change a given angle axis rotation p1 by the rotation p2
    """
    cdef Vec3 v3 = c_rotate_aa(to_vec(p1), to_vec(p2))
    return vec_to_np(v3)

cpdef mx2aa(mx):
    cdef Vec3 p = rot_mat_to_aa(to_mat(mx))
    return vec_to_np(p)

cpdef aa2q(aa):
    """convert angle axis to quaternion
    
    Parameters
    ----------
    aa: angle axis vector of lenth 3
    
    Returns
    -------
    Q: quaternion of length 4
    """
    return vec4_to_np(aa_to_quaternion(to_vec(aa)))
    