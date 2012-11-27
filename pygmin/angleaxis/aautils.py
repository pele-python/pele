import numpy as np
from pygmin.utils import rotations
from pygmin.angleaxis import CoordsAdapter
from pygmin.potentials.fortran.rmdrvt import rmdrvt as rotMatDeriv

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
        
class AASystem(object):
    ''' 
        Angle axis system wrapper
        
        Parameters
        ----------
            sites :  iteratable
                list of sites in the system
                
        Attributes
        ----------
        sites :
            list of all sites in the system
        
        
    '''
    def __init__(self, sites=[]):
        self.sites = sites
        
    def add_sites(self, sites):
        '''
            Add sites to the system
            
            Paramters
            ---------
            sites : iteratable
                list of sites
        '''
        self.sites.append(sites)
        
    def coords_adapter(self, coords=None):
        '''
        Create a coords adapter to easy access coordinate array
        
        Parameters
        ----------
        
        coords : ndarray, optional
            coords array to wrap
            
        '''
        return CoordsAdapter(nrigid=len(self.sites), coords=coords)
    
def aasitespring(com1, p1, com2, p2, sitetype=AASiteType()):
    '''
    calculate spring force between 2 sites
    
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
    returns: tuple
        spring cart, spring rot
    '''
    R1, R11, R12, R13 = rotMatDeriv(p1, True)
    R2 = rotations.aa2mx(p2)
    dR = R2 - R1
    
    dR = dR
    
    g_M = -2.*sitetype.M*(com2-com1)
    # dR_kl S_lm dR_km
    g_P = np.zeros(3) 
    g_P[0] = -2.*np.trace(np.dot(R11, np.dot(sitetype.S, dR.transpose()))) 
    g_P[1] = -2.*np.trace(np.dot(R12, np.dot(sitetype.S, dR.transpose()))) 
    g_P[2] = -2.*np.trace(np.dot(R13, np.dot(sitetype.S, dR.transpose()))) 

    return g_M, g_P

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

def aadistance(aasystem, coords1, coords2):
    ca1 = aasystem.coords_adapter(coords=coords1)
    ca2 = aasystem.coords_adapter(coords=coords2)
    
    d_sq = 0
    # first distance for sites only
    for i in xrange(ca1.nrigid):
        d_sq += aasitedistance_sq(ca1.posRigid[i], ca1.rotRigid[i],
                                  ca2.posRigid[i], ca2.rotRigid[i],
                                  aasystem.sites[i])
             
    return np.sqrt(d_sq)

def aainterpolate(initial, final, t):
    conf = initial.copy()
    for i in xrange(conf.shape[0]):
        conf[i] = rotations.q2aa(rotations.q_slerp(rotations.aa2q(initial[i]),
                                rotations.aa2q(final[i]), t))
        
# calculate the spring force on x1 to x2
def aaspringforce(aasystem, coords1, coords2):
    ca1 = aasystem.coords_adapter(coords=coords1)
    ca2 = aasystem.coords_adapter(coords=coords2)
    spring = aasystem.coords_adapter(np.zeros(coords1.shape))
    
    # first distance for sites only
    for i in xrange(ca1.nrigid):
        g_M, g_P = aasitespring(ca1.posRigid[i], ca1.rotRigid[i],
                                  ca2.posRigid[i], ca2.rotRigid[i],
                                  aasystem.sites[i])
        spring.posRigid[i] += g_M
        spring.rotRigid[i] += g_P
        
    return spring.coords

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
    print aasitespring(X1, p1, X2, p2, sitetype)
    g_M = np.zeros(3)
    g_P = np.zeros(3)
    
    for i in xrange(3):
        eps = 1e-6
        delta = np.zeros(3)
        delta[i] = eps
        g_M[i] = (aasitedistance_sq(X1+delta, p1, X2, p2, sitetype) 
                  - aasitedistance_sq(X1, p1, X2, p2, sitetype))/eps
        g_P[i] = (aasitedistance_sq(X1, p1+delta, X2, p2, sitetype) 
                  - aasitedistance_sq(X1, p1, X2, p2, sitetype))/eps
    print g_M, g_P
    xx = aasitespring(X1, p1, X2, p2, sitetype)
    print g_M/xx[0], g_P/xx[1]