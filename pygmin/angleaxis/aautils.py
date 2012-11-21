import numpy as np
from pygmin.utils import rotations
from pygmin.angleaxis import CoordsAdapter

__all__ = ["AASiteType", "aasitedistance_sq", "aadistance"]

class AASiteType(object):
    '''
    Definition of a angle axis site
    
    Parameters
    ----------
    
    M : float
        mass of the angle axis site
    I : 3x3 array
        momentum of inertia tensor
        
    '''  
    def __init__(self, M=1., S=np.identity(3, dtype="float64")):
        self.M = M
        self.S = S
        
        
def aasitedistance_sq(com1, p1, com2, p2, sitetype=AASiteType()):
    '''
    distance measure between 2 angle axis bodies of same type
    
    Parameters
    ----------
    com1:
        center of mass of 1st site
    p1: 
        angle axis vector of 1st site
    com2:
        center of mass of 2nd site
    p2:
        angle axis vector of 2nd site    
    sitetype: AASiteType, optional
        amgle axis site type with mass and moment of inertia tensor
    returns:
        distance squared
    '''
    R1 = rotations.aa2mx(p1)
    R2 = rotations.aa2mx(p2)
    dR = R2 - R1
    
    dR = dR
    
    d_M = sitetype.M*np.sum((com2-com1)**2)
    # dR_kl S_lm dR_km 
    d_P = np.trace(np.dot(dR, np.dot(sitetype.S, dR.transpose()))) 

    return d_M + d_P

def aadistance(aasites, coords1, coords2):
    ca1 = CoordsAdapter(coords=coords1)
    ca2 = CoordsAdapter(coords=coords2)
    
    d_sq = 0
    # first distance for sites only
    for i in xrange(ca1.nrigid):
        d_sq += aasitedistance_sq(ca1.posRigid[i], ca1.rotRigid[i],
                                  ca2.posRigid[i], ca2.rotRigid[i],
                                  aasites[i])     
    return np.sqrt(d_sq)

if __name__ == "__main__":
    natoms = 8
    x = np.random.random([natoms,3])*5
    x -= np.sum(x, axis=0)/natoms
    S=np.zeros([3,3])
    for i in xrange(3):
        for j in xrange(3):
            S[i][j]=np.sum(x[:,i]*x[:,j])
    sitetype = AASiteType(M=natoms, S=S)
    X1=np.zeros(3)
    p1=np.zeros(3)
    X1 = np.random.random(3)*2.
    X2 = np.random.random(3)*2.
    p1 = rotations.random_aa()
    p2 = rotations.random_aa()
    R1 = rotations.aa2mx(p1)
    R2 = rotations.aa2mx(p2)
    
    x1 = np.dot(R1, x.transpose()).transpose() + X1
    x2 = np.dot(R2, x.transpose()).transpose() + X2
    print "site representation:", np.sum((x1-x2)**2)
    print aasitedistance_sq(X1, p1, X2, p2, sitetype)
    