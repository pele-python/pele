"""
Functions related to rotations

Most of these were adapted from victor's rotations.f90.  Not all functions from
the file have been implemented.
Warning, they have not all been tested in this format.

.. currentmodule:: pele.utils.rotations

.. autosummary::
    :toctree: generated/
    
    
    q_multiply
    aa2q
    q2aa
    q2mx
    mx2q
    mx2aa
    aa2mx
    random_q
    random_aa
    takestep_aa
    rotate_aa
    small_random_aa
    vec_random
    vec_random_ndim
    vector_random_uniform_hypersphere
    q_slerp

"""
import numpy as np
from pele.utils._cpp_utils import rotate_aa, mx2aa, aa2q, aa2mx, \
    rot_mat_derivatives

rot_epsilon = 1e-6
dprand = np.random.rand


def q_multiply(q0, q1):
    """multiply 2 quaternions q1, q2"""
    q3 = np.zeros(4)
    q3[0] = q0[0] * q1[0] - q0[1] * q1[1] - q0[2] * q1[2] - q0[3] * q1[3]
    q3[1] = q0[0] * q1[1] + q0[1] * q1[0] + q0[2] * q1[3] - q0[3] * q1[2]
    q3[2] = q0[0] * q1[2] - q0[1] * q1[3] + q0[2] * q1[0] + q0[3] * q1[1]
    q3[3] = q0[0] * q1[3] + q0[1] * q1[2] - q0[2] * q1[1] + q0[3] * q1[0]
    return q3


def q2aa(qin):
    """
    quaternion to angle axis
    
    Parameters
    ----------
    Q: quaternion of length 4
    
    Returns
    -------
    output V: angle axis vector of lenth 3
    """
    q = np.copy(qin)
    if q[0] < 0.: q = -q
    if q[0] > 1.0: q /= np.sqrt(np.dot(q, q))
    theta = 2. * np.arccos(q[0])
    s = np.sqrt(1. - q[0] * q[0])
    if s < rot_epsilon:
        p = 2. * q[1:4]
    else:
        p = q[1:4] / s * theta
    return p


def q2mx(qin):
    """quaternion to rotation matrix"""
    Q = qin / np.linalg.norm(qin)
    RMX = np.zeros([3, 3], np.float64)
    Q2Q3 = Q[1] * Q[2]
    Q1Q4 = Q[0] * Q[3]
    Q2Q4 = Q[1] * Q[3]
    Q1Q3 = Q[0] * Q[2]
    Q3Q4 = Q[2] * Q[3]
    Q1Q2 = Q[0] * Q[1]

    RMX[0, 0] = 2. * (0.5 - Q[2] * Q[2] - Q[3] * Q[3])
    RMX[1, 1] = 2. * (0.5 - Q[1] * Q[1] - Q[3] * Q[3])
    RMX[2, 2] = 2. * (0.5 - Q[1] * Q[1] - Q[2] * Q[2])
    RMX[0, 1] = 2. * (Q2Q3 - Q1Q4)
    RMX[1, 0] = 2. * (Q2Q3 + Q1Q4)
    RMX[0, 2] = 2. * (Q2Q4 + Q1Q3)
    RMX[2, 0] = 2. * (Q2Q4 - Q1Q3)
    RMX[1, 2] = 2. * (Q3Q4 - Q1Q2)
    RMX[2, 1] = 2. * (Q3Q4 + Q1Q2)
    return RMX


def mx2q(mi):
    """convert a rotation matrix to a quaternion
    
    see discussion at 
    http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
    """
    q = np.zeros(4)
    m = np.transpose(mi)  # simply because I copied it from fortran code.
    trace = m[0, 0] + m[1, 1] + m[2, 2]

    if trace > 0.:
        s = np.sqrt(trace + 1.0) * 2.0
        q[0] = 0.25 * s
        q[1] = (m[1, 2] - m[2, 1]) / s
        q[2] = (m[2, 0] - m[0, 2]) / s
        q[3] = (m[0, 1] - m[1, 0]) / s
    elif (m[0, 0] > m[1, 1]) and (m[0, 0] > m[2, 2]):
        s = np.sqrt(1.0 + m[0, 0] - m[1, 1] - m[2, 2]) * 2.0
        q[0] = (m[1, 2] - m[2, 1]) / s
        q[1] = 0.25 * s
        q[2] = (m[1, 0] + m[0, 1]) / s
        q[3] = (m[2, 0] + m[0, 2]) / s
    elif m[1, 1] > m[2, 2]:
        s = np.sqrt(1.0 + m[1, 1] - m[0, 0] - m[2, 2]) * 2.0
        q[0] = (m[2, 0] - m[0, 2]) / s
        q[1] = (m[1, 0] + m[0, 1]) / s
        q[2] = 0.25 * s
        q[3] = (m[2, 1] + m[1, 2]) / s
    else:
        s = np.sqrt(1.0 + m[2, 2] - m[0, 0] - m[1, 1]) * 2.0
        q[0] = (m[0, 1] - m[1, 0]) / s
        q[1] = (m[2, 0] + m[0, 2]) / s
        q[2] = (m[2, 1] + m[1, 2]) / s
        q[3] = 0.25 * s

    if q[0] < 0:
        q = -q
    
    return q


# js850> this is commented because it's not documented and seems to be the same as q2mx()
# def rot_q2mx(qin):
# m = np.zeros([3,3], np.float64)
# 
# q = qin / np.linalg.norm(qin)
# 
# sq = q**2
# 
# m[0,0] = ( sq[1] - sq[2] - sq[3] + sq[0])
# m[1,1] = (-sq[1] + sq[2] - sq[3] + sq[0])
#     m[2,2] = (-sq[1] - sq[2] + sq[3] + sq[0])
# 
#     tmp0 = q[1]*q[2]
#     tmp1 = q[0]*q[3]
#     m[1,0] = 2.0 * (tmp0 + tmp1)
#     m[0,1] = 2.0 * (tmp0 - tmp1)
# 
#     tmp0 = q[1]*q[3]
#     tmp1 = q[2]*q[0]
#     m[2,0] = 2.0 * (tmp0 - tmp1)
#     m[0,2] = 2.0 * (tmp0 + tmp1)
#     tmp0 = q[2]*q[3]
#     tmp1 = q[0]*q[1]
#     m[2,1] = 2.0 * (tmp0 + tmp1)
#     m[1,2] = 2.0 * (tmp0 - tmp1)
# 
#     return m


def random_q():
    """
    uniform random rotation in angle axis formulation
    
    Notes
    -----
    input: 3 uniformly distributed random numbers
    uses the algorithm given in
    K. Shoemake, Uniform random rotations, Graphics Gems III, pages 124-132. Academic, New York, 1992.
    This first generates a random rotation in quaternion representation. We should substitute this by
    a direct angle axis generation, but be careful: the angle of rotation in angle axis representation
    is NOT uniformly distributed
    """
    from numpy import sqrt, sin, cos, pi

    u = np.random.uniform(0, 1, [3])
    q = np.zeros(4, np.float64)
    q[0] = sqrt(1. - u[0]) * sin(2. * pi * u[1])
    q[1] = sqrt(1. - u[0]) * cos(2. * pi * u[1])
    q[2] = sqrt(u[0]) * sin(2. * pi * u[2])
    q[3] = sqrt(u[0]) * cos(2. * pi * u[2])
    return q


def random_aa():
    """return a uniformly distributed random angle axis vector"""
    return q2aa(random_q())


def takestep_aa(p, maxtheta):
    """change an angle axis vector by a small rotation"""
    p[:] = rotate_aa(p, small_random_aa(maxtheta))


def small_random_aa(maxtheta):
    """generate a small random rotation"""
    # first choose a random unit vector
    p = vec_random()

    # linear for too small steps
    # this is not completely right but should be ok
    if maxtheta < rot_epsilon:
        p = p * dprand() * maxtheta
        return p

    s = 1. / (np.sin(0.5 * maxtheta) ** 2)
    # now choose the angle theta in range 0:step
    # with distribution sin(0.5*theta)**2
    u = dprand() * maxtheta
    while dprand() > s * np.sin(0.5 * u) ** 2:
        u = dprand() * maxtheta
    p = p * u
    return p


def vec_random():
    """uniform random unit vector"""
    p = np.zeros(3)
    u1 = dprand()
    u2 = dprand()
    z = 2 * u1 - 1.
    p[0] = np.sqrt(1 - z * z) * np.cos(2. * np.pi * u2)
    p[1] = np.sqrt(1 - z * z) * np.sin(2. * np.pi * u2)
    p[2] = z
    return p


def vec_random_ndim(n):
    """n-dimensional uniform random unit vector"""
    v = np.random.normal(size=n)
    v /= np.linalg.norm(v)
    return v


def vector_random_uniform_hypersphere(k):
    """return a vector sampled uniformly in a hypersphere of dimension k"""
    if k == 3:
        # this function is much faster than the general one
        u = vec_random()
    else:
        u = vec_random_ndim(k)
    # draw the magnitude of the vector from a power law density:
    # draws samples in [0, 1] from a power distribution with positive exponent k - 1.
    p = np.random.power(k)
    return p * u


def q_slerp(a, b, t):
    if t <= 0.:
        return a
    if t >= 1.:
        return b
    costheta = np.dot(a, b)

    c = b
    # if theta > 180., go other direction
    if costheta < 0.0:
        costheta = -costheta
        c = -c

    #linear interpolate close to zero
    if costheta > 1.0 - 1e-5:
        return t * b + (1 - t) * b

    theta = np.arccos(costheta)

    return (np.sin((1.0 - t) * theta) * a + np.sin(t * theta) * c) / np.sin(theta)


#
# only testing below here
#


def test_vector_random_uniform_hypersphere():  # pragma: no cover
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt

    nvec = 1000
    r = np.zeros([nvec, 3])
    for i in range(nvec):
        r[i, :] = vector_random_uniform_hypersphere(3)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(r[:, 0], r[:, 1], r[:, 2])
    plt.show()


if __name__ == "__main__":
    test_vector_random_uniform_hypersphere()

