"""
Functions related to rotations

Most of these were adapted from victor's rotations.f90.  Not all functions from
the file have been implemented.
Warning, they have not all been tested in this format.
So far I have tested 
q2mx
q2aa
"""
import numpy as np
import copy

rot_epsilon = 1e-6

def q_multiply(q1, q2):
    """ multiply 2 quaternions q1, q2 """
    q3 = copy.copy(q1)

    q3[0] = q0[0]*q1[0]-q0[1]*q1[1]-q0[2]*q1[2]-q0[3]*q1[3]
    q3[1] = q0[0]*q1[1]+q0[1]*q1[0]+q0[2]*q1[3]-q0[3]*q1[2]
    q3[2] = q0[0]*q1[2]-q0[1]*q1[3]+q0[2]*q1[0]+q0[3]*q1[1]
    q3[3] = q0[0]*q1[3]+q0[1]*q1[2]-q0[2]*q1[1]+q0[3]*q1[0]
    return q3

def aa2q( AA ):
    """
    convert angle axis to quaternion
    input V: angle axis vector of lenth 3
    output Q: quaternion of length 4
    """
    q = np.zeros(4, np.float64)

    thetah = 0.5 * np.linalg.norm( AA ) 
    q[0]  = np.cos( thetah )

    # do linear expansion for small epsilon
    if thetah < rot_epsilon:
        q[1:] = 0.5 * AA
    else:
        q[1:] = 0.5 * np.sin(thetah) * AA / thetah

    # make sure to have normal form
    if q[0] < 0.0: q = -q
    return q

def q2aa( qin ):
    """
    quaternion to angle axis
    input Q: quaternion of length 4
    output V: angle axis vector of lenth 3
    """
    q = copy.copy(qin)
    if q[0] < 0.: q = -q
    if q[0] > 1.0: q /= np.sqrt(np.dot(q,q))
    theta = 2. * np.arccos(q[0])
    s = np.sqrt(1.-q[0]*q[0])
    if s < rot_epsilon:
        p = 2. * q[1:4]
    else:
        p = q[1:4] / s * theta

    return p


def q2mx( qin ):
    """quaternion to rotation matrix"""
    Q = qin / np.linalg.norm(qin)
    RMX = np.zeros([3,3], np.float64)
    Q2Q3 = Q[1]*Q[2];
    Q1Q4 = Q[0]*Q[3];
    Q2Q4 = Q[1]*Q[3];
    Q1Q3 = Q[0]*Q[2];
    Q3Q4 = Q[2]*Q[3];
    Q1Q2 = Q[0]*Q[1];

    RMX[0,0] = 2.*(0.5 - Q[2]*Q[2] - Q[3]*Q[3]);
    RMX[1,1] = 2.*(0.5 - Q[1]*Q[1] - Q[3]*Q[3]);
    RMX[2,2] = 2.*(0.5 - Q[1]*Q[1] - Q[2]*Q[2]);
    RMX[0,1] = 2.*(Q2Q3 - Q1Q4);
    RMX[1,0] = 2.*(Q2Q3 + Q1Q4);
    RMX[0,2] = 2.*(Q2Q4 + Q1Q3);
    RMX[2,0] = 2.*(Q2Q4 - Q1Q3);
    RMX[1,2] = 2.*(Q3Q4 - Q1Q2);
    RMX[2,1] = 2.*(Q3Q4 + Q1Q2);
    return RMX

def rot_q2mx(qin):
    m = np.zeros([3,3], np.float64)

    q = qin / np.linalg.norm(qin)

    sq = q**2

    m[0,0] = ( sq[1] - sq[2] - sq[3] + sq[0])
    m[1,1] = (-sq[1] + sq[2] - sq[3] + sq[0])
    m[2,2] = (-sq[1] - sq[2] + sq[3] + sq[0])

    tmp0 = q[1]*q[2]
    tmp1 = q[0]*q[3]
    m[1,0] = 2.0 * (tmp0 + tmp1)
    m[0,1] = 2.0 * (tmp0 - tmp1)

    tmp0 = q[1]*q[3]
    tmp1 = q[2]*q[0]
    m[2,0] = 2.0 * (tmp0 - tmp1)
    m[0,2] = 2.0 * (tmp0 + tmp1)
    tmp0 = q[2]*q[3]
    tmp1 = q[0]*q[1]
    m[2,1] = 2.0 * (tmp0 + tmp1)
    m[1,2] = 2.0 * (tmp0 - tmp1)

    return m

def aa2mx( p ):
    return q2mx( aa2q( p ) )
