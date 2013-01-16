import numpy as np
from pygmin.utils import rotations
from pygmin.angleaxis import CoordsAdapter
from pygmin.potentials.fortran.rmdrvt import rmdrvt as rotMatDeriv
from pygmin.transition_states import interpolate_linear
from math import pi
from pygmin import takestep
from pygmin.transition_states import zeroev
from pygmin.angleaxis.aamindist import TransformAngleAxisCluster

__all__ = ["AASiteType", "AASystem", "interpolate_angleaxis"]

def interpolate_angleaxis(initial, final, t):
    ''' interpolate between 2 arrays of angle axis coordinates
    
    Parameters
    ----------    
    initial: np.array
        configuration 1
    final: np.array
        configuration 2
    t: float
        interpolation parameter with t in [0,1]    
    
    '''    
    conf = initial.copy()
    for i in xrange(conf.shape[0]):
        conf[i] = rotations.q2aa(rotations.q_slerp(rotations.aa2q(initial[i]),
                                rotations.aa2q(final[i]), t))
    return conf


class AASiteType(object):
    '''
    Definition of a angle axis site
    
    Parameters
    ----------
    
    M : float
        total mass of the angle axis site
    S : 3x3 array
        weighted tensor of gyration S_ij = \sum m_i x_i x_j 
        
    '''
    def __init__(self, M=1., S=np.identity(3, dtype="float64"), cog=np.zeros(3), W = 1.0):
        self.M = M
        self.S = S
        self.cog = cog
        self.W = W
        
        self.inversion = None
        self.symmetries = []
            
    def distance_squared(self, com1, p1, com2, p2):
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
        
        d_M = self.W*np.sum((com2-com1)**2)
        # dR_kl S_lm dR_km 
        d_P = np.trace(np.dot(dR, np.dot(self.S, dR.transpose()))) 
        d_mix = 2.*self.W * np.dot((com2-com1), np.dot(dR, self.cog))
        return d_M + d_P + d_mix

    def distance_squared_grad(self, com1, p1, com2, p2):
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
        
        g_M = -2.*self.W*(com2-com1)
        # dR_kl S_lm dR_km
        g_P = np.zeros(3) 
        g_P[0] = -2.*np.trace(np.dot(R11, np.dot(self.S, dR.transpose())))
        g_P[1] = -2.*np.trace(np.dot(R12, np.dot(self.S, dR.transpose())))
        g_P[2] = -2.*np.trace(np.dot(R13, np.dot(self.S, dR.transpose())))
    
        g_M -= 2.*self.W *  np.dot(dR, self.cog)
        g_P[0] -= 2.*self.W * np.dot((com2-com1), np.dot(R11, self.cog))
        g_P[1] -= 2.*self.W * np.dot((com2-com1), np.dot(R12, self.cog))
        g_P[2] -= 2.*self.W * np.dot((com2-com1), np.dot(R13, self.cog))

        return g_M, g_P

        
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
                list of AASiteType
        '''
        self.sites += sites
        
    def coords_adapter(self, coords=None):
        '''
        Create a coords adapter to easy access coordinate array
        
        Parameters
        ----------
        
        coords : ndarray, optional
            coords array to wrap
            
        '''
        return CoordsAdapter(nrigid=len(self.sites), coords=coords)
    

    def distance_squared(self, coords1, coords2):
        ''' Calculate the squared distance between 2 configurations'''
        ca1 = self.coords_adapter(coords=coords1)
        ca2 = self.coords_adapter(coords=coords2)
        
        d_sq = 0
        # first distance for sites only
        for i in xrange(ca1.nrigid):
            d_sq += self.sites[i].distance_squared(
                                           ca1.posRigid[i], ca1.rotRigid[i],
                                           ca2.posRigid[i], ca2.rotRigid[i])
                 
        return np.sqrt(d_sq)
                
    # calculate the spring force on x1 to x2
    def distance_squared_grad(self, coords1, coords2):
        ''' Calculate gradient with respect to coords 1 for the squared distance'''
        ca1 = self.coords_adapter(coords=coords1)
        ca2 = self.coords_adapter(coords=coords2)
        spring = self.coords_adapter(np.zeros(coords1.shape))
        
        # first distance for sites only
        for i in xrange(ca1.nrigid):
            g_M, g_P = self.sites[i].distance_squared_grad(
                                           ca1.posRigid[i], ca1.rotRigid[i],
                                           ca2.posRigid[i], ca2.rotRigid[i])
            spring.posRigid[i] += g_M
            spring.rotRigid[i] += g_P
            
        return spring.coords
    
    def neb_distance(self, coords1, coords2, distance=True, grad=True):
        d=None
        if distance:
            d = self.distance_squared(coords1, coords2)
        g = None
        if grad:
            g = self.distance_squared_grad(coords1, coords2)
        return d, g
    
    def interpolate(self, initial, final, t):
        cinitial = self.coords_adapter(initial)
        cfinal = self.coords_adapter(final)
        
        cnew = self.coords_adapter(np.zeros_like(initial))
        cnew.coords[:] = interpolate_linear(initial, final, t)
        cnew.posRigid[:] = interpolate_linear(cinitial.posRigid, cfinal.posRigid, t)
        cnew.rotRigid[:] = interpolate_angleaxis(cinitial.rotRigid, cfinal.rotRigid, t)        
        return cnew.coords
  
    def align_coords(self, x1, x2):
        c2 = self.coords_adapter(x1)
        c1 = self.coords_adapter(x2)
        for p1, p2 in zip(c1.rotRigid,c2.rotRigid):
            if np.linalg.norm(p2) < 1e-6:
                if(np.linalg.norm(p1) < 1e-6):
                    continue
                n2 = p1/np.linalg.norm(p1)*2.*pi
            else:
                n2 = p2/np.linalg.norm(p2)*2.*pi
        
            while True:
                p2n = p2+n2
                if(np.linalg.norm(p2n - p1) > np.linalg.norm(p2 - p1)):
                    break
                p2[:]=p2n
                
            while True:
                p2n = p2-n2
                if(np.linalg.norm(p2n - p1) > np.linalg.norm(p2 - p1)):
                    break
                p2[:]=p2n 
    
    def align_path(self, path):
        for i in xrange(1, len(path)):
            c2 = self.coords_adapter(path[i])
            c1 = self.coords_adapter(path[i-1])
            for p1, p2 in zip(c1.rotRigid,c2.rotRigid):
                if np.linalg.norm(p2) < 1e-6:
                    if(np.linalg.norm(p1) < 1e-6):
                        continue
                    n2 = p1/np.linalg.norm(p1)*2.*pi
                else:
                    n2 = p2/np.linalg.norm(p2)*2.*pi
            
                while True:
                    p2n = p2+n2
                    if(np.linalg.norm(p2n - p1) > np.linalg.norm(p2 - p1)):
                        break
                    p2[:]=p2n
                    
                while True:
                    p2n = p2-n2
                    if(np.linalg.norm(p2n - p1) > np.linalg.norm(p2 - p1)):
                        break
                    p2[:]=p2n
                    
    def zeroEV(self, x):
        zev = []
        ca = self.coords_adapter(x)
        translate = zeroev.zeroEV_translation(ca.posRigid)
        #rotate_r = zeroev.zeroEV_rotation(ca.posRigid)
        #rotate_aa = 
        transform = TransformAngleAxisCluster(self)
        d = 1e-6
        dx = x.copy()
        transform.rotate(dx, rotations.aa2mx(np.array([d, 0, 0])))
        dx -= x
        dx /= d
        
        dy = x.copy()
        transform.rotate(dy, rotations.aa2mx(np.array([0, d, 0])))
        dy -= x
        dy /= d
        
        dz = x.copy()
        transform.rotate(dz, rotations.aa2mx(np.array([0, 0, d])))
        dz -= x
        dz /= d
        
        return zev + [dx, dy, dz]
    
    def orthogopt(self, v, coords):
        zev = self.zeroEV(coords)
        zeroev.orthogonalize(v, zev)
        return v
    
    def ortogopt_aa(self, v, coords):
        v = v.copy()
        zev = zeroev.gramm_schmidt(self.zeroEV(coords))
        zeroev.orthogonalize(v, zev)
        return v

class TakestepAA(takestep.TakestepInterface):
    def __init__(self, topology, rotate=1.6, translate=0.):
        self.rotate = rotate
        self.translate = translate
        self.topology = topology
        self
    
    def takeStep(self, coords, **kwargs):
        ca = self.topology.coords_adapter(coords)
        takestep.uniform_displace(self.translate, ca.posRigid)
        takestep.rotate(self.rotate, ca.rotRigid)

    def scale(self, factor):
        self.translate *= factor
        self.rotate *= factor
        
        
if __name__ == "__main__":
    natoms = 3
    x = np.random.random([natoms,3])*5
    masses = [1., 1., 16.] #np.random.random(natoms)
    print masses
    x -= np.average(x, axis=0, weights=masses)
    cog = np.average(x, axis=0)
    S=np.zeros([3,3])
    for i in xrange(3):
        for j in xrange(3):
            S[i][j]=np.sum(x[:,i]*x[:,j])
    site = AASiteType(M=natoms, S=S, W=natoms, cog=cog)
    X1=np.zeros(3)
    p1=np.zeros(3)    
    
    X1 = 10.1*np.random.random(3)
    X2 = 10.1*np.random.random(3)
    p1 = rotations.random_aa()
    p2 = rotations.random_aa()
        
    R1 = rotations.aa2mx(p1)
    R2 = rotations.aa2mx(p2)
    
    x1 = np.dot(R1, x.transpose()).transpose() + X1
    x2 = np.dot(R2, x.transpose()).transpose() + X2
    
    print "site representation:", np.sum((x1-x2)**2)
    print "distance function:  ", site.distance_squared(X1, p1, X2, p2)
    print site.distance_squared_grad(X1, p1, X2, p2)
    g_M = np.zeros(3)
    g_P = np.zeros(3)
    
    for i in xrange(3):
        eps = 1e-6
        delta = np.zeros(3)
        delta[i] = eps
        g_M[i] = (site.distance_squared(X1+delta, p1, X2, p2) 
                  - site.distance_squared(X1, p1, X2, p2))/eps
        g_P[i] = (site.distance_squared(X1, p1+delta, X2, p2) 
                  - site.distance_squared(X1, p1, X2, p2))/eps
    print g_M, g_P
    xx = site.distance_squared_grad(X1, p1, X2, p2)
    print g_M/xx[0], g_P/xx[1]
