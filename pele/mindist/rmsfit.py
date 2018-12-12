from __future__ import print_function
import numpy as np
from pele.utils import rotations

__all__ = ["findrotation", "findrotation_kabsch", "findrotation_kearsley"]

def findrotation_kabsch(coords1, coords2, align_com=True):
    """
    Kabsch, Wolfgang, (1976) "A solution of the best rotation to relate two sets of vectors", Acta Crystallographica 32:922
    
    ..note::
        this has different return values than findrotation_kearsley.  The return values for this
        function may change in the future.
    """
    # check if arrays are of same size
    if coords1.size != coords2.size:
        raise ValueError("dimension of arrays does not match")
    
    # reshape the arrays
    x1 = coords1.reshape([-1,3]).copy()
    x2 = coords2.reshape([-1,3]).copy()
    
    # determine number of atoms
    natoms = x1.shape[0]
    
    # set both com to zero
    if align_com:
        com1 = np.sum(x1, axis=0) / float(natoms)
        com2 = np.sum(x2, axis=0) / float(natoms)
        x1 -= com1
        x2 -= com2
  
    # calculate covariance matrix
    A = np.dot( x2.transpose(), x1)
    # and do single value decomposition
    u, s, v = np.linalg.svd(A)
 
    if np.linalg.det(u) * np.linalg.det(v) + 1.0 < 1e-8:
        s[-1] = -s[-1]
        u[:,-1] = -u[:,-1]

    return  np.dot(u, v).transpose()
    
def findrotation_kearsley(x1, x2, align_com=True):
    """Return the rotation matrix which aligns XB with XA
    
    Return the matrix which
    aligns structure XB to be as similar as possible to structure XA.
    To be precise, rotate XB, so as to minimize the distance |XA - XB|.

    Rotations will be done around the origin, not the center of mass

    Rotational alignment follows the prescription of
    Kearsley, Acta Cryst. A, 45, 208-210, 1989
    http://dx.doi.org/10.1107/S0108767388010128
    """
    if x1.size != x2.size:
        raise ValueError("dimension of arrays does not match")
    
    # reshape the arrays
    x1 = x1.reshape([-1,3]).copy()
    x2 = x2.reshape([-1,3]).copy()
    # determine number of atoms
    natoms = x1.shape[0]
    
    # set both com to zero
    if align_com:
        com1 = np.sum(x1,axis=0) / float(natoms)
        com2 = np.sum(x2,axis=0) / float(natoms)
        x1 -= com1
        x2 -= com2

    x1 = x1.ravel() 
    x2 = x2.ravel()
    
    # TODO: this is very dirty!
    #########################################
    # Create matrix QMAT
    #########################################

    QMAT = np.zeros([4,4], np.float64)
    for J1 in range(natoms):
        J2 = 3* J1 -1
        XM = x1[J2+1] - x2[J2+1]
        YM = x1[J2+2] - x2[J2+2]
        ZM = x1[J2+3] - x2[J2+3]
        XP = x1[J2+1] + x2[J2+1]
        YP = x1[J2+2] + x2[J2+2]
        ZP = x1[J2+3] + x2[J2+3]
        QMAT[0,0] = QMAT[0,0] + XM**2 + YM**2 + ZM**2
        QMAT[0,1] = QMAT[0,1] - YP*ZM + YM*ZP
        QMAT[0,2] = QMAT[0,2] - XM*ZP + XP*ZM
        QMAT[0,3] = QMAT[0,3] - XP*YM + XM*YP
        QMAT[1,1] = QMAT[1,1] + YP**2 + ZP**2 + XM**2
        QMAT[1,2] = QMAT[1,2] + XM*YM - XP*YP
        QMAT[1,3] = QMAT[1,3] + XM*ZM - XP*ZP
        QMAT[2,2] = QMAT[2,2] + XP**2 + ZP**2 + YM**2
        QMAT[2,3] = QMAT[2,3] + YM*ZM - YP*ZP
        QMAT[3,3] = QMAT[3,3] + XP**2 + YP**2 + ZM**2

    QMAT[1,0] = QMAT[0,1]
    QMAT[2,0] = QMAT[0,2]
    QMAT[2,1] = QMAT[1,2]
    QMAT[3,0] = QMAT[0,3]
    QMAT[3,1] = QMAT[1,3]
    QMAT[3,2] = QMAT[2,3]

    ###########################################
    """
    Find eigenvalues and eigenvectors of QMAT.  The eigenvector corresponding
    to the smallest eigenvalue is the quaternion which rotates XB into best
    alignment with XA.  The smallest eigenvalue is the squared distance between
    the resulting structures.
    """
    ###########################################
    (eigs, vecs) = np.linalg.eig(QMAT)

    imin = np.argmin(eigs)
    eigmin = eigs[imin] # the minimum eigenvector
    Q2 = vecs[:,imin]  # the eigenvector corresponding to the minimum eigenvalue
    if eigmin < 0.:
        if abs(eigmin) < 1e-6:
            eigmin = 0.
        else:
            print('minDist> WARNING minimum eigenvalue is ',eigmin,' change to absolute value')
            eigmin = -eigmin

    dist = np.sqrt(eigmin) # this is the minimized distance between the two structures

    Q2 = np.real_if_close(Q2, 1e-10)
    if np.iscomplexobj(Q2):
        raise ValueError("Q2 is complex")
    return dist, rotations.q2mx(Q2)

findrotation = findrotation_kearsley

if __name__ == "__main__":
    x1 = np.random.random(24)
    mx = rotations.q2mx(rotations.random_q())
    
    x2 = np.dot(mx,x1.reshape(-1,3).transpose()).transpose().reshape(-1)
    print(x2-x1)
    print(mx-findrotation_kabsch(x1,x2))
    print(findrotation_kabsch(x1,x2))
    print(findrotation_kearsley(x1,x2)[1])
    

