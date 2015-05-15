"""
# distutils: language = C++

This module provides access to the c++ routines that deal with rotations
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
    Vec3 c_rotate_aa "pele::rotate_aa" (Vec3 &, Vec3 &) except +
    Vec3 rot_mat_to_aa(Matrix33 & mx) except +
    Vec4 aa_to_quaternion(Vec3 & mx) except +
    Matrix33 aa_to_rot_mat(Vec3 & mx) except +
    
    void c_rot_mat_derivatives "pele::rot_mat_derivatives"(
            Vec3 & p,
            Matrix33 & rmat,
            Matrix33 & drm1,
            Matrix33 & drm2,
            Matrix33 & drm3) except +

cdef Vec3 to_vec(p) except *:
    cdef Vec3 v
    p = np.asarray(p)
    if p.size != 3:
        raise ValueError("p must be a 3 element array or list")
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

cdef Vec4 to_vec4(p) except *:
    cdef Vec4 v
    p = np.asarray(p)
    if p.size != 4:
        raise ValueError("p must be a 4 element array or list")
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
    A = np.asarray(A).ravel()
    if A.size != 9:
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
    return Anew.reshape([3,3])

cpdef rotate_aa(p1, p2):
    """change a given angle axis rotation p2 by the rotation p1
    """
    return vec_to_np(c_rotate_aa(to_vec(p1), to_vec(p2)))

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

cpdef aa2mx(aa):
    """convert an angle axis rotation to a rotation matrix"""
    return mat_to_np(aa_to_rot_mat(to_vec(aa)))

cpdef rot_mat_derivatives(p, with_grad=True):
    """compute the derivatives of a rotation matrix
    
    Parameters
    ----------
    p : ndarray
        angle axis vector
    
    Returns
    -------
    R : the rotation matrix corresponding to p
    dR1, dR2, dR3 : derivatives of R in the x, y, z directions.
    """
    if not with_grad:
        return aa2mx(p)
    cdef Matrix33 R, dR1, dR2, dR3
    c_rot_mat_derivatives(to_vec(p), R, dR1, dR2, dR3)
    return mat_to_np(R), mat_to_np(dR1), mat_to_np(dR2), mat_to_np(dR3)
    
    
