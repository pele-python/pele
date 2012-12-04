import numpy as np
import aautils

class RigidFragment(aautils.AASiteType):
    ''' defines a single rigid fragment 
    
    In the most simple case, this is just a whole molecule
    '''
    
    def __init__(self):
        aautils.AASiteType.__init__(self)
        self.atom_positions = []
        self.atom_types = []
        self.atom_masses = []
        
    def add_atom(self, atomtype, pos, mass=1.0):
        '''Add a new atom to the rigid fragment
        
        Parameters
        ----------
        type: string
            type identifier
        pos: np.array
            position of the atom
        mass: mass of the atom
        '''
        self.atom_types.append(atomtype)
        self.atom_positions.append(pos)
        self.atom_masses.append(mass)
        
    def finalize_setup(self):
        '''finalize setup after all sites have been added
        
        This will shift the center of mass to the origin and calculate
        the total mass and weighted tensor of gyration
        '''
        
        # first correct for the center of mass
        com = np.average(self.atom_positions, axis=0, weights=self.atom_masses)
        for x in self.atom_positions:
            x[:] -= com
        
        # calculate total mass
        self.M = np.sum(self.atom_masses)
        
        # now calculate the weighted moment of inertia tensor
        self.S[:] = 0.
        for x, m in zip(self.atom_positions, self.atom_masses):
            self.S[:] += m*np.outer(x, x)
            
    def to_atomistic(self, com, R):
        return com + np.dot(R, np.transpose(self.atom_positions)).transpose()
            
    def transform_grad(self, p, g):
        g_com = np.sum(g, axis=0)
        R, R1, R2, R3 = aautils.rotMatDeriv(p, True)
        g_p = np.zeros_like(g_com)
        g_p[0] = np.sum(
                        np.dot(g, np.dot(R1, np.transpose(self.atom_positions))),
                        )
        g_p[1] = np.sum(
                        np.dot(g, np.dot(R2, np.transpose(self.atom_positions)).transpose()),
                        )
        g_p[2] = np.sum(
                        np.dot(g, np.dot(R3, np.transpose(self.atom_positions)).transpose()),
                        )
        return g_com, g_p

    def redistribute_forces(self, p, grad_com, grad_p):
        R, R1, R2, R3 = aautils.rotMatDeriv(p, True)
        grad = np.dot(R1, np.transpose(self.atom_positions)).transpose()*grad_p[0]
        grad += np.dot(R2, np.transpose(self.atom_positions)).transpose()*grad_p[1]
        grad += np.dot(R3, np.transpose(self.atom_positions)).transpose()*grad_p[2]
        grad += grad_com
        
        grad = (self.atom_masses * grad.transpose()).transpose()/self.M
        
        return grad

            
class RBSystem(aautils.AASystem):
    def __init__(self):
        aautils.AASystem.__init__(self)
        self.indices=[]
        self.natoms=0
        
    def add_sites(self, sites):
        aautils.AASystem.add_sites(self, sites)
        
        for site in sites:
            nsite_atoms = len(site.atom_positions)
            if not hasattr(site, "atom_indices"):
                site.indices = range(self.natoms, self.natoms+nsite_atoms)
            self.natoms += nsite_atoms
        
    def to_atomistic(self, rbcoords):
        ca = self.coords_adapter(rbcoords)
        atomistic = np.zeros([self.natoms,3])
        for site, com, p in zip(self.sites, ca.posRigid, ca.rotRigid):
            atoms = site.to_atomistic(com, p)
            for i,x in zip(site.indices, atoms):
                atomistic[i]=x
        return atomistic

    def transform_grad(self, rbcoords, grad):
        ca = self.coords_adapter(rbcoords)
        rbgrad = self.coords_adapter(np.zeros_like(rbcoords))
        
        for site, p, g_com, g_p in zip(self.sites, ca.rotRigid,
                                       rbgrad.posRigid,rbgrad.rotRigid):
            g_com[:], g_p[:] = site.transform_grad(p, grad[site.indices])
        return rbgrad.coords
                
    def redistribute_forces(self, rbcoords, rbgrad):
        ca = self.coords_adapter(rbcoords)
        cg = self.coords_adapter(rbgrad)
        grad = np.zeros([self.natoms,3])
        for site, p, g_com, g_p in zip(self.sites, ca.rotRigid, cg.posRigid, cg.rotRigid):
            gatom = site.redistribute_forces(p, g_com, g_p)
            for i,x in zip(site.indices, gatom):
                grad[i]=x
        return grad
    
if __name__ == "__main__":
    from math import sin, cos, pi
    from copy import deepcopy
    water = RigidFragment()
    rho   = 0.9572
    theta = 104.52/180.0*pi      
    water.add_atom("O", np.array([0., 0., 0.]), 16.)
    water.add_atom("H", rho*np.array([0.0, sin(0.5*theta), cos(0.5*theta)]), 1.)
    water.add_atom("H", rho*np.array([0.0, -sin(0.5*theta), cos(0.5*theta)]), 1.)
    water.finalize_setup()
    # define the whole water system
    system = RBSystem()
    nrigid = 1
    system.add_sites([deepcopy(water) for i in xrange(nrigid)])

    rbcoords = np.random.random(6*nrigid)
    print rbcoords
    coords = system.to_atomistic(rbcoords)
    print coords
    grad = np.random.random(coords.shape)
    rbgrad = system.transform_grad(rbcoords, grad)
    print rbgrad
    gnew = system.redistribute_forces(rbcoords, rbgrad)
    print gnew
    print system.transform_grad(rbcoords, gnew)