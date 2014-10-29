import numpy as np

from pele.angleaxis import aatopology
from pele.potentials.potential import potential
from pele.mindist import ExactMatchAtomicCluster
from pele.utils.rotations import mx2aa
from pele.utils import rotations

class RigidFragment(aatopology.AASiteType):
    '''Defines a single rigid fragment 
    
    Notes
    -----
    This defines a collection of atoms that compose a single rigid bodies.
    This class collects all the information necessary to perform operations
    on the rigid body such as converting from center of mass + angle axis rotataion
    to atomistic coordinates. In the most simple case, this is just a whole molecule
    '''
    
    def __init__(self):
        aatopology.AASiteType.__init__(self)
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
            position of the atom  sn402: I think this is in absolute Cartesian coords?
        mass: mass of the atom
        '''
        self.atom_types.append(atomtype)
        self.atom_positions.append(pos.copy())
        self.atom_masses.append(mass)
        
    def finalize_setup(self, shift_com=True):
        '''finalize setup after all sites have been added
        
        This will shift the center of mass to the origin and calculate
        the total mass and weighted tensor of gyration
        '''
        # sn402: This doesn't need to be altered for PBCs, but bear in mind that 
        # self.atom_positions can be negative

        # first correct for the center of mass
        com = np.average(self.atom_positions, axis=0, weights=self.atom_masses)
        if shift_com:
            for x in self.atom_positions:
                x[:] -= com

        self.cog = np.average(self.atom_positions, axis=0)
        self.W = float(len(self.atom_masses))
        
        # calculate total mass
        self.M = np.sum(self.atom_masses)
        
        # now calculate the weighted moment of inertia tensor
        self.S[:] = 0.
        for x in self.atom_positions:
            self.S[:] += np.outer(x, x)
        self.Sm = np.zeros([3,3])
        for x, m in zip(self.atom_positions, self.atom_masses):
            self.Sm[:] += m*np.outer(x, x)
        self._determine_symmetries()  # sn402: check this.
        
        # calculate aa rotations for later
        self.symmetriesaa = []
        for rot in self.symmetries:
            self.symmetriesaa.append(mx2aa(rot))
            
    def to_atomistic(self, com, p):
        """convert the center of mass position + angle axis vector to atomistic coords
        """
        R = rotations.aa2mx(p)
        return com + np.dot(R, np.transpose(self.atom_positions)).transpose()
            
    def transform_grad(self, p, g):
        g_com = np.sum(g, axis=0)
        R, R1, R2, R3 = rotations.rot_mat_derivatives(p)
        g_p = np.zeros_like(g_com)
        for ga, x in zip(g, self.atom_positions):
            g_p[0] += np.dot(ga, np.dot(R1, x))
            g_p[1] += np.dot(ga, np.dot(R2, x))
            g_p[2] += np.dot(ga, np.dot(R3, x))
#                        

#        g_p[0] = -np.sum(
#                        np.dot(g, np.dot(R1, np.transpose(self.atom_positions)).transpose())
#                        )
#        g_p[1] = np.sum(
#                        np.dot(g, np.dot(R2, np.transpose(self.atom_positions)).transpose())
#                        )
#        g_p[2] = np.sum(
#                        np.dot(g, np.dot(R3, np.transpose(self.atom_positions)).transpose())
#                        )
        return g_com, g_p

    def redistribute_forces(self, p, grad_com, grad_p):
        R, R1, R2, R3 = rotations.rot_mat_derivatives(p)
        grad = np.dot(R1, np.transpose(self.atom_positions)).transpose()*grad_p[0]
        grad += np.dot(R2, np.transpose(self.atom_positions)).transpose()*grad_p[1]
        grad += np.dot(R3, np.transpose(self.atom_positions)).transpose()*grad_p[2]
        grad += grad_com
        
        grad = (self.atom_masses * grad.transpose()).transpose()/self.M
        
        return grad


    def _determine_inversion(self, permlist):
        x = np.array(self.atom_positions).flatten()
        xi = -x
        exact = ExactMatchAtomicCluster(permlist=permlist, can_invert=False)
        
        for rot, invert in exact.standard_alignments(x, xi):
            if exact.check_match(x, xi, rot, invert):
                self.inversion = exact._last_checked_rotation
                return 

    def _determine_rotational_symmetry(self, permlist):
        x = np.array(self.atom_positions).flatten()
        exact = ExactMatchAtomicCluster(permlist=permlist, can_invert=False)
        for rot, invert in exact.standard_alignments(x, x):
            if exact.check_match(x, x, rot, invert):
                rot = exact._last_checked_rotation
                exists=False
                for rot2 in self.symmetries:
                    if np.linalg.norm(rot2 - rot) < 1e-6:
                        exists=True
                        break
                if not exists:
                    self.symmetries.append(rot)
                
        
    def _determine_symmetries(self):
        self.symmetries = []
        self.inversion = None
        
        perm_dict = dict()
        for t in self.atom_types:
            perm_dict[t] = []
            
        for i, t in zip(xrange(len(self.atom_types)), self.atom_types):
            perm_dict[t].append(i)
            
        permlist = []
        for i in perm_dict.itervalues():
            if len(i) > 1:
                permlist.append(i)
        
        self._determine_inversion(permlist)
        self._determine_rotational_symmetry(permlist)
        
# sn402 Added this.        
class RigidFragmentBulk(RigidFragment):
    """Modified site type class; this is an exact copy of RigidFragment 
       but with an overloaded method "get_smallest_rij" to incorporate 
       periodic boundary conditions.
    """
    def __init__(self,boxl):
        aatopology.AASiteType.__init__(self)
        self.atom_positions = []
        self.atom_types = []
        self.atom_masses = []
        self.boxsize = boxl
    
    def get_smallest_rij(self, com1, com2):
        """return the shortest vector from com1 to com2 (both numpy arrays containing 
        coordinates for any number of atoms) using periodic boundary conditions.
        """
        #print "Using periodic boundary conditions, Python version"
        boxvec = self.boxsize
        dx = com2 - com1  
        dx = dx.reshape(-1, boxvec.size) # sn402: takes dx (however long it is - 
        # in this case 3*nrigid) and reshapes it to boxvec.size columns.
        dx -= boxvec * np.round(dx / boxvec[np.newaxis,:]) # np.newaxis inserts a new column into boxvec. 
        # This is just to match shape with dx and hence allow division.
        return dx       
                    
class RBTopology(aatopology.AATopology):
    """This defines the topology of a collection of rigid bodies.
    """
    def __init__(self):
        aatopology.AATopology.__init__(self)
        self.natoms=0
        
    def get_atomtypes(self):
        atom_types = [None for i in xrange(self.natoms)]
        for site in self.sites:
            for i, atype in zip(site.atom_indices, site.atom_types):
                atom_types[i]=atype
        return atom_types

    def get_atommasses(self):
        masses = [None for i in xrange(self.natoms)]
        for site in self.sites:
            for i, mass in zip(site.atom_indices, site.atom_masses):
                masses[i]=mass
        return masses
        
    def add_sites(self, sites):
        """
        add sites to the Rigid Body topology
        
        Parameters
        ----------
        sites : RigidFragment
            usually this will be an object of type RigidFragment
        """
        aatopology.AATopology.add_sites(self, sites)
        for site in sites:
            nsite_atoms = len(site.atom_positions)
            if hasattr(site, "atom_indices"):
                print "warning: the c++ RBPotentialWrapper does not support user defined atom_indices.  The potential may be wrong"
            else:
                site.atom_indices = range(self.natoms, self.natoms+nsite_atoms)
            self.natoms += nsite_atoms
    
    def finalize_setup(self):
        from pele.angleaxis import _cpp_aa
        self.set_cpp_topology(_cpp_aa.cdefRBTopology(self))

            
    def get_atom_labels(self):
        labels=[]
        for s in self.sites:
            for t in s.atom_types:
                labels.append(str(t))
        return labels
    
    def to_atomistic(self, rbcoords):
        """convert rigid body coords to atomistic coords

        Note: there is a c++ implementation of which is much faster
        """
        ca = self.coords_adapter(rbcoords)
        atomistic = np.zeros([self.natoms,3])
        for site, com, p in zip(self.sites, ca.posRigid, ca.rotRigid):
            atoms = site.to_atomistic(com, p)
            for i,x in zip(site.atom_indices, atoms):
                atomistic[i]=x
        return atomistic

    def transform_gradient(self, rbcoords, grad):
        """convert atomistic gradient into a gradient in rigid body coords

        Note: there is a c++ implementation of which is much faster
        """

        ca = self.coords_adapter(rbcoords)
        rbgrad = self.coords_adapter(np.zeros_like(rbcoords))
        for site, p, g_com, g_p in zip(self.sites, ca.rotRigid,
                                       rbgrad.posRigid,rbgrad.rotRigid):
            g_com[:], g_p[:] = site.transform_grad(p, grad.reshape(-1,3)[site.atom_indices])
        return rbgrad.coords
                
    def redistribute_gradient(self, rbcoords, rbgrad):
        ca = self.coords_adapter(rbcoords)
        cg = self.coords_adapter(rbgrad)
        grad = np.zeros([self.natoms,3])
        for site, p, g_com, g_p in zip(self.sites, ca.rotRigid, cg.posRigid, cg.rotRigid):
            gatom = site.redistribute_forces(p, g_com, g_p)
            for i,x in zip(site.atom_indices, gatom):
                grad[i]=x
        return grad
    
        
    def set_cpp_topology(self, cpp_topology):
        """provide class to access the fast c++ topology routines"""
        self.cpp_topology = cpp_topology


# sn402: Added to transfer the PBC information
class RBTopologyBulk(RBTopology, aatopology.AATopologyBulk):
    def __init__(self, boxvec, sites=None):
        if sites is None:
            sites = []
        self.sites = sites
        self.boxvec = boxvec
        self.natoms=0
        
    
class RBPotentialWrapper(potential):
    """Wrap a potential
    
    Note: there is a c++ implementation of which is much faster
    """
    def __init__(self, rbsystem, pot):
        self.pot = pot
        self.rbsystem = rbsystem
        
    def getEnergy(self, rbcoords):
        coords = self.rbsystem.to_atomistic(rbcoords)
        return self.pot.getEnergy(coords.reshape(-1))
    
    def getEnergyGradient(self, rbcoords):
        coords = self.rbsystem.to_atomistic(rbcoords)
        E, g = self.pot.getEnergyGradient(coords.reshape(-1))
        return E, self.rbsystem.transform_gradient(rbcoords, g)   # This uses the metric tensor
    

def test(): # pragma: no cover
    from math import sin, cos, pi
    from copy import deepcopy
    water = RigidFragment()
    rho   = 0.9572
    theta = 104.52/180.0*pi      
    water.add_atom("O", np.array([0., -1., 0.]), 1.)
    water.add_atom("O", np.array([0., 1., 0.]), 1.)
    #water.add_atom("H", rho*np.array([0.0, sin(0.5*theta), cos(0.5*theta)]), 1.)
    #water.add_atom("H", rho*np.array([0.0, -sin(0.5*theta), cos(0.5*theta)]), 1.)
    water.finalize_setup()
    # define the whole water system
    system = RBTopology()
    nrigid = 1
    system.add_sites([deepcopy(water) for i in xrange(nrigid)])
    from pele.utils import rotations
    rbcoords=np.zeros(6)
    rbcoords[3:]= rotations.random_aa()
    
    coords = system.to_atomistic(rbcoords)
    
    print "rb coords\n", rbcoords
    print "coords\n", coords
    grad = (np.random.random(coords.shape) -0.5)
    
    v = coords[1] - coords[0]
    v/=np.linalg.norm(v)
    
    grad[0] -= np.dot(grad[0], v)*v
    grad[1] -= np.dot(grad[1], v)*v
    
    grad -= np.average(grad, axis=0)
    grad /= np.linalg.norm(grad)
    
    print "torque", np.linalg.norm(np.cross(grad, v))
    rbgrad = system.transform_gradient(rbcoords, grad)
    p = rbcoords[3:]
    x = rbcoords[0:3]
    gp = rbgrad[3:]
    gx = rbgrad[:3]
    
    R, R1, R2, R3 = rotations.rot_mat_derivatives(p)        
    
    print "test1", np.linalg.norm(R1*gp[0])     
    print "test2", np.linalg.norm(R2*gp[1])     
    print "test3", np.linalg.norm(R3*gp[2])
    print "test4", np.linalg.norm(R1*gp[0]) + np.linalg.norm(R2*gp[1]) + np.linalg.norm(R3*gp[2])
         
    dR = R1*gp[0] + R2*gp[1] + R3*gp[2]
    print "test", np.linalg.norm(R1*gp[0] + R2*gp[1] + R3*gp[2])
    print np.trace(np.dot(dR, dR.transpose()))
    #G = water.metric_tensor_aa(p)
    #print np.dot(p, np.dot(G, p))     
    exit()
    
    
    gnew = system.redistribute_forces(rbcoords, rbgrad)
    print gnew
    print system.transform_grad(rbcoords, gnew)
    
def test_bulk_class():
    
    boxvec = np.array([5,10,20])
    coords1 = np.array([1,2,3,4,4,4])
    coords2 = np.array([4,8,8,-4,-4,-4])
    print boxvec, coords1, coords2
    a = RBTopologyBulk(boxvec)
    for i in range(2):
        otp = RigidFragmentBulk(boxvec)   # sn402: changed
        otp.add_atom("O", np.array([0.0, 0.0, 0.0]), 1.)
        otp.add_atom("O", np.array([1.0, 0.0, 0.0]), 1.)        
        otp.finalize_setup()     
        a.add_sites(otp)
   
    b = a.distance_squared(coords1, coords2)
    print b      
    
if __name__ == "__main__":
    #test()
    test_bulk_class()