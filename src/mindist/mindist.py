import numpy as np
import rotations as rot

def alignCoM( X1, X2):
    """
    align the center of mass of X2 with that of X1
    """
    natoms = len(X1) / 3
    for i in xrange(3):
        com = np.sum( X1[i::3] - X2[i::3] )
        com /= natoms
        X2[i::3] += com

def CoMToOrigin( X1):
    """
    move the center of mass to the origin
    """
    natoms = len(X1) / 3
    for i in xrange(3):
        com = np.sum( X1[i::3] )
        com /= natoms
        X1[i::3] -= com
    return X1


def getAlignRotation(XA, XB):
    """
    Return the quaternion which
    aligns structure XB to be as similar as possible to structure XA.
    To be precise, rotate XB, so as to minimize the distance |XA - XB|.

    Rotations will be done around the origin, not the center of mass

    Rotational alignment follows the prescription of
    Kearsley, Acta Cryst. A, 45, 208-210, 1989
    http://dx.doi.org/10.1107/S0108767388010128
    """
    nsites = len(XA)/3

    #move the center of mass to the origin
    #CoMToOrigin(XA)
    #CoMToOrigin(XB)

    #distcm = np.linalg.norm(XA-XB)
    #print "distcm", distcm

    #########################################
    #Create matrix QMAT
    #########################################

    QMAT = np.zeros([4,4], np.float64)
    for J1 in xrange(nsites):
       J2 = 3*(J1) -1
       XM = XA[J2+1] - XB[J2+1]
       YM = XA[J2+2] - XB[J2+2]
       ZM = XA[J2+3] - XB[J2+3]
       XP = XA[J2+1] + XB[J2+1]
       YP = XA[J2+2] + XB[J2+2]
       ZP = XA[J2+3] + XB[J2+3]
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
    #print "eigenvalues", eigs

    imin = np.argmin(eigs)
    eigmin = eigs[imin] #the minimum eigenvector
    Q2 = vecs[:,imin]  #the eigenvector corresponding to the minimum eigenvalue
    if eigmin < 0.:
        if abs(eigmin) < 1e-6:
            eigmin = 0.
        else:
            print 'minDist> WARNING minimum eigenvalue is ',eigmin,' change to absolute value'
            eigmin = -eigmin

    dist = np.sqrt(eigmin) #this is the minimized distance between the two structures
    #print "dist from eigenvalue", dist
    #print "Q2", Q2, "norm", np.linalg.norm(Q2)
    #aa = rot.q2aa( Q2)
    #print "aa ", aa, "norm", np.linalg.norm(aa)

    return dist, Q2

def alignRotation(XA, XB):
    """
    Align structure XB to be as similar as possible to structure XA.
    To be precise, rotate XB, so as to minimize the distance |XA - XB|.

    Rotations will be done around the origin, not the center of mass
    """
    nsites = len(XA)/3

    dist, Q2 = getAlignRotation(XA, XB)
    ###################################################################
    # Q2 contains the quaternion which rotates XB to best align with X1.
    # rotate XB according to Q2
    ###################################################################

    rot_mx = rot.q2mx( Q2 )
    #print rot_mx
    for j in range(nsites):
        i = 3*j
        XB[i:i+3] = np.dot( rot_mx, XB[i:i+3] )

    return dist


def minDist(X1, X2):
    """
    Minimize the distance between two clusters.  The following symmetries will be accounted for
    
    Translational symmetry

    Global rotational symmetry
    """
    #alignCoM(X1, X2)
    X1 = CoMToOrigin(X1)
    X2 = CoMToOrigin(X2)

    #align rotation degrees of freedom
    dist = alignRotation(X1, X2)
    return dist


def main():
    natoms = 5
    X1 = np.random.uniform(-1,1,[natoms*3])*(float(natoms))**(1./3)
    X2 = np.random.uniform(-1,1,[natoms*3])*(float(natoms))**(1./3)

    #X1 = np.array( [ 0., 0., 0., 1., 0., 0., 0., 0., 1.,] )
    #X2 = np.array( [ 0., 0., 0., 1., 0., 0., 0., 1., 0.,] )
    import copy
    X1i = copy.copy(X1)
    X2i = copy.copy(X2)

    distinit = np.linalg.norm(X1-X2)
    print "distinit", distinit

    dist = minDist(X1,X2)
    distfinal = np.linalg.norm(X1-X2)
    print "dist from eigenvalue", dist
    print "distfinal", distfinal

    import printing.print_atoms_xyz as printxyz
    with open("out.xyz", "w") as fout:
        CoMToOrigin(X1i)
        CoMToOrigin(X2i)
        printxyz.printAtomsXYZ(fout, X1i )
        printxyz.printAtomsXYZ(fout, X2i )
        printxyz.printAtomsXYZ(fout, X1 )
        printxyz.printAtomsXYZ(fout, X2 )

if __name__ == "__main__":
    main()
