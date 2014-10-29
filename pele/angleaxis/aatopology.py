"""this module holds the classes applicable for working
with general rigid body systems.  i.e. those that do not
necessarily have a representation as a set of atomistic coords.
see rigidbody.py for those classes which derive from these.
"""

import numpy as np
from pele.utils import rotations
from pele.angleaxis import CoordsAdapter
from pele.transition_states import interpolate_linear
from math import pi
from pele import takestep
from pele.transition_states import _zeroev as zeroev
from pele.angleaxis.aamindist import TransformAngleAxisCluster

from pele.utils.rotations import rot_mat_derivatives
from _aadist import sitedist_grad, sitedist

__all__ = ["AASiteType", "AATopology", "interpolate_angleaxis", "TakestepAA"]

def interpolate_angleaxis(initial, final, t):
    '''interpolate between 2 arrays of angle axis coordinates
    
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
    '''Definition of an angle axis site
    
    Each angle axis site is fully characterized by a tensor of gyration S for the shape,
    and the center of geometry, which can differ from the center of mass.
    TODO: see paper to be published
    
    Parameters
    ----------
    
    M : float
        total mass of the angle axis site
    W : float
        sum of all weights
    S : 3x3 array
        weighted tensor of gyration S_ij = \sum m_i x_i x_j 
    cog : 3 dim np.array
        center of gravity
    inversion : 3x3 np.array
        matrix which defines on how to perform an inversion which doesn't
        affect the coordinates. None if no inversion can be performed
    symmetries : list of 3x3 np.arrays
        list of all symmetry operations which can be performed on the angle axis site
        excluding an inversion
        
    '''
    def __init__(self, M=1., S=np.identity(3, dtype="float64"), cog=np.zeros(3), W = 1.0, Inv=None, Sym=None):
        self.M = M
        self.S = S
        self.cog = cog
        self.W = W
        
        self.inversion = None
        if Sym is None:
            self.symmetries = [np.eye( 3 ), ]

    def get_smallest_rij(self, com1, com2):
        """return the shortest vector from com1 to com2 (both numpy arrays containing 
        coordinates for any number of atoms)
        """
        return com2 - com1
          
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
            angle axis site type with mass and moment of inertia tensor
        returns:
            distance squared
        '''
        print self.get_smallest_rij(com1, com2)
        return sitedist(self.get_smallest_rij(com1, com2), p1, p2, self.S, self.W, self.cog)

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
        
        return sitedist_grad(self.get_smallest_rij(com1, com2), p1, p2, self.S, self.W, self.cog)
    
    def metric_tensor(self, p):
        '''calculate the mass weighted metric tensor '''
        R, R1, R2, R3 = rot_mat_derivatives(p, True)        
        g = np.zeros([3,3])
                
        g[0,0] = np.trace(np.dot(R1, np.dot(self.Sm, R1.transpose())))
        g[0,1] = np.trace(np.dot(R1, np.dot(self.Sm, R2.transpose())))
        g[0,2] = np.trace(np.dot(R1, np.dot(self.Sm, R3.transpose())))
        g[1,1] = np.trace(np.dot(R2, np.dot(self.Sm, R2.transpose())))
        g[1,2] = np.trace(np.dot(R2, np.dot(self.Sm, R3.transpose())))
        g[2,2] = np.trace(np.dot(R3, np.dot(self.Sm, R3.transpose())))
        
        g[1,0] = g[0,1]
        g[2,1] = g[1,2]
        g[2,0] = g[0,2]
        gx = np.identity(3)*self.M        

        return gx, g
        
    def metric_tensor_cog(self, x, p):
        '''calculate the metric tensor when for w_i != m_i '''
        R, R1, R2, R3 = rot_mat_derivatives(p, True)        
        g = np.zeros([6,6])

        # the com part
        g[0:3,0:3] = self.W * np.identity(3)

        # the rotational part        
        g[3,3] = np.trace(np.dot(R1, np.dot(self.Sm, R1.transpose())))
        g[3,4] = np.trace(np.dot(R1, np.dot(self.Sm, R2.transpose())))
        g[3,5] = np.trace(np.dot(R1, np.dot(self.Sm, R3.transpose())))
        g[4,4] = np.trace(np.dot(R2, np.dot(self.Sm, R2.transpose())))
        g[4,5] = np.trace(np.dot(R2, np.dot(self.Sm, R3.transpose())))
        g[5,5] = np.trace(np.dot(R3, np.dot(self.Sm, R3.transpose())))
        
        g[4,3] = g[3,4]
        g[5,4] = g[4,5]
        g[5,3] = g[3,5]
        
        # the mixing part
#        g[0:3,3] = 2.*self.W * np.dot(R1, self.cog)
#        g[0:3,4] = 2.*self.W * np.dot(R2, self.cog)
#        g[0:3,5] = 2.*self.W * np.dot(R3, self.cog)
#        #g[0:3,3:] = g[0:3,3:].transpose() 
#        
#        g[3,0:3] = g[0:3,3] 
#        g[4,0:3] = g[0:3,4] 
#        g[5,0:3] = g[0:3,5] 
        
        return g
    
class AATopology(object):
    ''' 
        Angle axis topology
        
        An angle axis topology stores all topology information for an angle axis 
        system. The AATopology is composed of several angle axis sites (AASiteType),
        which describe the shape of the angle axis site and each site carries a 
        position and orientation. Therefore, the length of the coordinate array
        must be 6*number_of_sites.
        
        Parameters
        ----------
            sites :  iteratable
                list of sites in the system
                
        Attributes
        ----------
        sites :
            list of all sites in the system
        
        
    '''
    def __init__(self, sites=None):
        if sites is None:
            sites = []
        self.sites = sites
        
    def add_sites(self, sites):
        '''
            Add a site to the topology
            
            Parameters
            ---------
            sites : iteratable
                list of AASiteType
        '''
        self.sites += sites
        
    def coords_adapter(self, coords=None):
        ''' Create a coords adapter to easy access coordinate array '''
        return CoordsAdapter(nrigid=len(self.sites), coords=coords)
    

    def _distance_squared_python(self, coords1, coords2):
        ''' Calculate the squared distance between 2 configurations'''

        ca1 = self.coords_adapter(coords=coords1)
        ca2 = self.coords_adapter(coords=coords2)
        
        d_sq = 0
        # first distance for sites only
        for i in xrange(ca1.nrigid):
            d_sq += self.sites[i].distance_squared(
                                           ca1.posRigid[i], ca1.rotRigid[i],
                                           ca2.posRigid[i], ca2.rotRigid[i])
# sn402: change the distance function between sites
# sn402: looks like it is actually ok
        return d_sq
    
    def distance_squared(self, coords1, coords2):
        '''Calculate the squared distance between 2 configurations'''
        try:
            return self.cpp_topology.distance_squared(coords1, coords2)
        except AttributeError:
            return self._distance_squared_python(coords1, coords2)
    
     
    def _distance_squared_grad_python(self, coords1, coords2):
        '''Calculate gradient with respect to coords 1 for the squared distance'''
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
    
    def distance_squared_grad(self, coords1, coords2):
        '''Calculate gradient with respect to coords 1 for the squared distance'''
        try:
            return self.cpp_topology.distance_squared_grad(coords1, coords2)
        except AttributeError:
            return self._distance_squared_grad_python(coords1, coords2)
    
    def neb_distance(self, coords1, coords2, distance=True, grad=True):
        '''wrapper function called by neb to get distance between 2 images '''
        d = None
        if distance:
            d = self.distance_squared(coords1, coords2)
        g = None
        if grad:
            g = self.distance_squared_grad(coords1, coords2)
        return d, g
    
    def interpolate(self, initial, final, t):
        ''' 
        interpolate between 2 sets of angle axis configurations 
        
        initial : np.array
            initial configuration
        final : np.array
            final configuration
        t : float
            interpolation parameters [0,1]
        '''
        cinitial = self.coords_adapter(initial)
        cfinal = self.coords_adapter(final)
        
        cnew = self.coords_adapter(np.zeros_like(initial))
        cnew.coords[:] = interpolate_linear(initial, final, t)
        cnew.posRigid[:] = interpolate_linear(cinitial.posRigid, cfinal.posRigid, t)
        cnew.rotRigid[:] = interpolate_angleaxis(cinitial.rotRigid, cfinal.rotRigid, t)        
        return cnew.coords

# js850> 9/2014 I commented this out because it seemed like it wasn't used 
#               and the documentation suggests it shouldn't be used
#    def align_coords(self, x1, x2):
#        ''' align the angle axis coordinates to minimize |p2 - p1|
#        
#            The angle axis vectors are perodic, this function changes the definition
#            of the angle axis vectors in x2 to match closest the ones in x1. This can
#            be useful if simple distances are used on angle axis vectors. However,
#            using the difference in angle axis vectors should be avoided, instead
#            the angle axis distance function should be used which properly takes care
#            of the rotatational degrees of freedoms
#        '''
#        c2 = self.coords_adapter(x1)
#        c1 = self.coords_adapter(x2)
#        for p1, p2 in zip(c1.rotRigid,c2.rotRigid):
#            if np.linalg.norm(p2) < 1e-6:
#                if(np.linalg.norm(p1) < 1e-6):
#                    continue
#                n2 = p1/np.linalg.norm(p1)*2.*pi
#            else:
#                n2 = p2/np.linalg.norm(p2)*2.*pi
#        
#            while True:
#                p2n = p2+n2
#                if(np.linalg.norm(p2n - p1) > np.linalg.norm(p2 - p1)):
#                    break
#                p2[:]=p2n
#                
#            while True:
#                p2n = p2-n2
#                if(np.linalg.norm(p2n - p1) > np.linalg.norm(p2 - p1)):
#                    break
#                p2[:]=p2n 
    
    def align_path(self, path):
        """ensure a series of images are aligned with each other
        
        Parameters
        ----------
        path : list of arrays
            This is a list of numpy array in com + angle axis format
        
        Notes
        -----
        this simply aligns the angle axis vectors
        """
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
    
    def _zeroEV_python(self, x):
        """return a list of zero eigenvectors
        
        This does both translational and rotational eigenvectors
        """
        zev = []
        ca = self.coords_adapter(x)
        cv = self.coords_adapter(np.zeros(x.shape))
            
        # get the zero eigenvectors corresponding to translation
        translate_rigid = zeroev.zeroEV_translation(ca.posRigid)
        
        for v in translate_rigid:
            cv.posRigid[:] = v
            zev.append(cv.coords.copy())
        
        # get the zero eigenvectors corresponding to rotation
        #rotate_r = zeroev.zeroEV_rotation(ca.posRigid)
        #rotate_aa = 
        transform = TransformAngleAxisCluster(self)
        d = 1e-5
        dx = x.copy()
        transform.rotate(dx, rotations.aa2mx(np.array([d, 0, 0])))
        self.align_path([x, dx])
        dx -= x
        dx /= np.linalg.norm(dx)
        
        dy = x.copy()
        transform.rotate(dy, rotations.aa2mx(np.array([0, d, 0])))
        self.align_path([x, dy])
        dy -= x
        dy /= np.linalg.norm(dy)
        
        dz = x.copy()
        transform.rotate(dz, rotations.aa2mx(np.array([0, 0, d])))
        self.align_path([x, dz])
        dz -= x
        dz /= np.linalg.norm(dz)
        
        #print "Zero eigenvectors", zev         
        return zev + [dx, dy, dz]
    
    def zeroEV(self, x):
        """return a list of zero eigenvectors
        
        This does both translational and rotational eigenvectors
        """
        try:
            return self.cpp_topology.get_zero_modes(x)
        except AttributeError:
            return self._zeroEV_python(x)
    
    def orthogopt(self, v, coords):
        """make v orthogonal to the zero eigenvectors at position coords""" 
        zev = zeroev.gramm_schmidt(self.zeroEV(coords))
        zeroev.orthogonalize(v, zev)
        return v
    
# js850> 9/2014 I commented this out because it is not used, not documented, and spelled wrong
#    def ortogopt_aa(self, v, coords):
#        v = v.copy()
#        zev = zeroev.gramm_schmidt(self.zeroEV(coords))
#        zeroev.orthogonalize(v, zev)
#        return v
    
    def metric_tensor(self, coords):
        '''get the metric tensor for a current configuration '''
        ca = self.coords_adapter(coords=coords)
        g = np.zeros([coords.size, coords.size])
        offset = 3*ca.nrigid
        # first distance for sites only
        for i in xrange(ca.nrigid):
            g_M, g_P = self.sites[i].metric_tensor(ca.rotRigid[i])
            g[3*i:3*i+3, 3*i:3*i+3] = g_M
            g[3*i+offset:3*i+3+offset, 3*i+offset:3*i+3+offset] = g_P
            
        return g
          

# sn402: new class to call the correct (PBC) versions of the cpp distance functions.          
class AATopologyBulk(AATopology):
    
    def __init__(self, boxvec, sites=None):
        if sites is None:
            sites = []
        self.sites = sites
        self.boxvec = boxvec  #sn402: should probably support the case where no boxvec is passed in
        
    def distance_squared(self, coords1, coords2):
        '''Calculate the squared distance between 2 configurations'''
        try:   
            return self.cpp_topology.distance_squared_bulk(coords1, coords2, self.boxvec)
        except AttributeError:
            print "Warning: used Python version of AATopologyBulk.distance_squared"
            return self._distance_squared_python(coords1, coords2)

            
    def distance_squared_grad(self, coords1, coords2):
        '''Calculate gradient with respect to coords 1 for the squared distance'''
        try:
            return self.cpp_topology.distance_squared_grad_bulk(coords1, coords2, self.boxvec)
        except AttributeError:
            print "Warning: used Python version of AATopologyBulk.distance_squared_grad"            
            return self._distance_squared_grad_python(coords1, coords2)
        
          
class TakestepAA(takestep.TakestepInterface):
    def __init__(self, topology, rotate=1.6, translate=0.):
        self.rotate = rotate
        self.translate = translate
        self.topology = topology
    
    def takeStep(self, coords, **kwargs):
        ca = self.topology.coords_adapter(coords)
        takestep.uniform_displace(self.translate, ca.posRigid)
        takestep.rotate(self.rotate, ca.rotRigid)

    def scale(self, factor):
        self.translate *= factor
        self.rotate *= factor
        
        
def test(): # pragma: no cover
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
    
    import _aadist
    print "site representation:", np.sum((x1-x2)**2)
    print "distance function:  ", site.distance_squared(X1, p1, X2, p2)

    print "fortran function:  ", _aadist.sitedist(X2 - X1, p1, p2, site.S, site.W, cog)

    coords1 = np.random.random(120)
    coords2 = np.random.random(120)
    
    import time
    t0 = time.time()
    for i in xrange(1000):
        site.distance_squared(X1, p1, X2, p2)
    t1 = time.time()
    print "time python", t1-t0
    for i in xrange(1000):
        sitedist(X2 - X1, p1, p2, site.S, site.W, cog)
    
 #_aadist.aadist(coords1, coords2, site.S, site.W, cog)
    t2 = time.time()
    print "time fortran", t2-t1
#    for i in xrange(1000/20):
#        #_aadist.sitedist(X1, p1, X2, p2, site.S, site.W, cog)
#        _aadist.aadist(coords1, coords2, site.S, site.W, cog)
    t2 = time.time()
    print "time fortran acc", t2-t1
    
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
    print _aadist.sitedist_grad(X2 - X1, p1, p2, site.S, site.W, cog)
#    print _aadist.sitedist_grad(com1, p1, com2, p2, self.S, self.W, self.cog)



if __name__ == "__main__":
    test()
