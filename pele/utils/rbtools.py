"""
.. currentmodule:: pele.utils.rbtools
.. autosummary:: 
    :toctree: generated/

    pele.utils.rbtools.CoordsAdapter 

"""
from numpy import vectorize

__all__ = ["CoordsAdapter"]
class CoordsAdapter(object):
    '''Wrapper to access coordinate array for rigid body systems
    
    The CoordsAdapter is a wrapper for a coords array. It creates views to directly
    access rigid body position & rotations, atom positions and lattice coordinates.
    This offers a convenient way to access coorinates without the hazzle of 
    indexing

    :Example:

    >>> import numpy as np
    >>> from pele.utils.rbtools import CoordsAdapter
    >>> 
    >>> nrigid = 10
    >>> coords = np.zeros(6*nrigid)
    >>> ca = CoordsAdapter(nrigid=nrigid, coords=coords)
    >>>
    >>> ca.posRigid[0] = np.random.random(3)
    >>> ca.rotRigid[0] = np.random.random(3)
    
    '''
    
    nrigid = 0
    ''' number of rigid bodies '''
    natoms = 0
    ''' number of single atoms '''
    nlattice = 0
    ''' number of lattice degrees of freedom '''
    coords = None
    ''' coordinate array '''
    
    posAtoms=None
    ''' array view for atom positions of dimension [3,nrigid] '''
    posRigid=None
    ''' array view for rigid body positions of dimension [3,nrigid] '''
    rotRigid=None
    ''' array view for rigid body rotations of dimension [3,nrigid] '''
    lattice=None
    ''' array view for lattice coordinates of dimension [nlattice] '''
    
    def __init__(self, nrigid=None, natoms=None, nlattice=0, coords=None):
        ''' initialise the coordinate wrapper
        
        Initialises the coordinate wrapper. The coordinates array can be
        either specified directly in the constructor or later changed
        via updateCoords.
        
        :param nrigid: number of rigid bodies
        :type nrigid: int
        :param natoms: number of single atoms
        :type natoms: int
        :param nlattice: number of lattice degrees of freedom
        :type nlattice: 0
        :param coords: the coordinate array
        :type coords: numpy.array
        '''
        if nrigid is None and natoms is None:
            nrigid = coords.size/6
            natoms = 0
            
        self.nrigid = nrigid
        self.natoms = natoms
        self.nlattice = nlattice
        if coords is not None:
            self.updateCoords(coords)

    def copy(self):
        return CoordsAdapter(nrigid=self.nrigid, natoms=self.natoms, nlattice=self.nlattice, coords=self.coords)
    
    def updateCoords(self, coords):
        ''' update the coordinate array
        
        This function can be called if the coordinate adapter should point
        to a different coordinate array.
        
        :param coords: the coordinate array
        :type coords: numpy.array
        '''
        natoms = self.natoms
        nrigid = self.nrigid
        self.coords = coords
        
        self.posAtoms=None
        self.posRigid=None
        self.rotRigid=None
        self.lattice=None
        
        if natoms > 0:
            self.posAtoms = self.coords[6*nrigid:6*nrigid+3*natoms].reshape(natoms, 3)
        
        if nrigid > 0:
            self.posRigid = self.coords[0:3*nrigid].reshape(nrigid, 3)
            self.rotRigid = self.coords[3*nrigid:6*nrigid].reshape(nrigid, 3)
    
        if self.nlattice > 0:
            self.lattice = self.coords[-self.nlattice:]


def test_com():
    import numpy as np
    from math import pi, sin, cos, atan2
    y = 2.1
    coords = np.array([-2.,-y,0.,-2,y,0,2.,-y,0.,2.,y,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.])
    boxvec = np.array([8,8,2])
    
    ca = CoordsAdapter(coords=coords)
    ndim = np.size(boxvec)
    print ca.nrigid
    print ca.posRigid    
    
    theta = 2.*pi*ca.posRigid/boxvec[np.newaxis,:]
    print theta        
    make_xi = np.vectorize(lambda x: cos(x))
    make_zeta = np.vectorize(lambda x: sin(x))            
    xi = make_xi(theta)
    zeta = make_zeta(theta)   
    print xi 
    print zeta         
    xi_ave = xi.sum(0)/ca.nrigid
    zeta_ave = zeta.sum(0)/ca.nrigid
    theta_ave = np.zeros(ndim)
    for i in range(ndim):
        theta_ave[i] = atan2(-zeta_ave[i],-xi_ave[i]) + pi
    print xi_ave 
    print zeta_ave 
    print theta_ave
    com = (theta_ave*boxvec/(2.*pi))%boxvec
    print com
    
if __name__ == "__main__":
    test_com()