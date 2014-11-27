from itertools import izip

import numpy as np

from pele.angleaxis import RBTopology, RBSystem, RigidFragment, RBPotentialWrapper
from pele.potentials import BasePotential
from pele.utils import rotations
from plate_potential import PlatePotential

EDGE1_TYPE = "O" 
EDGE2_TYPE = "C"
EDGE3_TYPE = "N"
OTHER_TYPE = "H"

def draw_pymol(coords):
    import pele.utils.pymolwrapper as pym
    pym.start()
    pym.draw_spheres(coords, "A", 1)

#class HarmonicPotential(BasePotential):
#    def __init__(self, atoms1, atoms2):
#        self.atoms1 = np.array(atoms1)
#        self.atoms2 = np.array(atoms2)
#    
#    def getEnergy(self, x):
#        e, g = self.getEnergyGradient(x)
#        return e
#    
#    def getEnergyGradient(self, x):
#        x = x.reshape([-1,3])
#        grad = np.zeros(x.shape)
#        etot = 0.
#        
#        for a1, a2 in izip(self.atoms1, self.atoms2):
#            dx = x[a1,:] - x[a2,:]
#            r2 = np.sum(dx**2)
#            etot += 0.5 * r2
#            grad[a1,:] += dx 
#            grad[a2,:] -= dx
#        return etot, grad.ravel()
            
        
class MolAtomIndexParser(object):
    """this tool helps with getting the correct indices for the atoms on a given edge of a plate"""
    def __init__(self, aatopology, nrigid):
        self.aatopology = aatopology
        atomtypes = self.aatopology.get_atomtypes()
        self.atom_types = np.array(atomtypes).reshape([nrigid, -1])
        self.atoms_per_mol = self.atom_types.shape[1]
    
    def get_atom_indices(self, molecule_number, atom_type):
        mol_atom_types = self.atom_types[molecule_number, :]
        atoms = np.where(mol_atom_types == atom_type)[0]
        
        atoms += molecule_number * self.atoms_per_mol
        atoms.sort()
        return list(atoms)
        
        
        

#class CombinePotential(BasePotential):
#    def __init__(self, potentials):
#        self.potentials = potentials
#    
#    def getEnergy(self, coords):
#        e = 0
#        for pot in self.potentials:
#            e += pot.getEnergy(coords)
#        return e
#    
#    def getEnergyGradient(self, coords):
#        etot = 0
#        grad = np.zeros(coords.size)
#        for pot in self.potentials:
#            e, g = pot.getEnergyGradient(coords)
#            etot += e
#            grad += g
#        return etot, grad.ravel()

def make_triangular_plate(atoms_per_side=8):
    """construct a single triangular plate
    """
    theta = 60. * np.pi / 180.
    v1 = np.array([1,0,0])
    v2 = np.array([0.5, np.sin(theta), 0])
    
    aps = atoms_per_side
    plate = RigidFragment()

    for i in xrange(aps-1):
        for j in xrange(aps-1):
            if i + j >= aps-1:
                break
            xnew = v1*i + v2*j
            if (i == 0 and j == 0 or 
                i == 0 and j == aps-2 or
                i == aps-2 and j == 0):
                atomtype = OTHER_TYPE
            elif i == 0:
                atomtype = EDGE1_TYPE
            elif j == 0:
                atomtype = EDGE2_TYPE
            elif i + j == aps-2:
                atomtype = EDGE3_TYPE
            else:
                atomtype = OTHER_TYPE
            plate.add_atom(atomtype, xnew, 1)

#    draw(coords)
    plate.finalize_setup()
    return plate


class PlateFolder(RBSystem):
    """
    This will build a system class for a cluster of interacting plates
    
    """
    def __init__(self, nmol):
        self.nrigid = nmol
        super(PlateFolder, self).__init__()
        
        self.setup_params(self.params)
        self._create_potential()

    def get_random_configuration(self):
        # js850> this is a bit sketchy because nrigid might not be defined here.
        # probably we can get the number of molecules some other way.
        coords = 10.*np.random.random(6*self.nrigid)
        ca = self.aasystem.coords_adapter(coords)
        for p in ca.rotRigid:
            p[:] = rotations.random_aa()
        return coords


    def setup_aatopology(self):
        """this sets up the topology for the whole rigid body system"""
        topology = RBTopology()
        topology.add_sites([make_triangular_plate() for i in xrange(self.nrigid)])
        
        self.render_scale = 0.2
        self.atom_types = topology.get_atomtypes()
        
        self.draw_bonds = []
#         for i in xrange(self.nrigid):
#             self.draw_bonds.append((3*i, 3*i+1))
#             self.draw_bonds.append((3*i, 3*i+2))
        
        topology.finalize_setup()
    
        return topology

    def setup_params(self, params):
        """set some system dependent parameters to imrprove algorithm performance"""
        
        params.double_ended_connect.local_connect_params.tsSearchParams.iprint = 10
        nebparams = params.double_ended_connect.local_connect_params.NEBparams
        nebparams.max_images = 50
        nebparams.image_density = 5
        nebparams.iter_density = 10.
        nebparams.k = 5.
        nebparams.reinterpolate = 50
        
        nebparams.NEBquenchParams["iprint"] = 10
        
        
        tssearch = params.double_ended_connect.local_connect_params.tsSearchParams
        tssearch.nsteps_tangent1 = 10
        tssearch.nsteps_tangent2 = 30
        tssearch.lowestEigenvectorQuenchParams["nsteps"] = 50
        tssearch.iprint=1
        tssearch.nfail_max = 100
        
        params.takestep.translate = 5.
    
    def _create_potential(self):
        """construct the potential which will compute the energy and gradient in atomistic (cartesian) coordinates
        
        The bonded (hinged) sides interact with an attractive harmonic potential.  Each atom 
        in the bond has a single interaction partner.
        
        The loosly attractive sides interact with an LJ potential.  These interactions are
        not specific.  Each atom interacts with every other one.
        
        All atoms repel each other with a WCA potential. 
        """
        parser = MolAtomIndexParser(self.aatopology, self.nrigid)
        
        # this is currently only set up for a tetrahedron
        assert self.nrigid == 4
        # do hinges
        harmonic_atoms1 = []
        harmonic_atoms2 = []
        harmonic_atoms1 += parser.get_atom_indices(0, EDGE1_TYPE)
        harmonic_atoms2 += parser.get_atom_indices(1, EDGE1_TYPE)
         
        harmonic_atoms1 += parser.get_atom_indices(0, EDGE2_TYPE)
        harmonic_atoms2 += parser.get_atom_indices(2, EDGE1_TYPE)
        
        harmonic_atoms1 += parser.get_atom_indices(0, EDGE3_TYPE)
        harmonic_atoms2 += parser.get_atom_indices(3, EDGE1_TYPE)
        
        self.harmonic_atoms = np.array(harmonic_atoms1 + harmonic_atoms2, np.integer)

        harmonic_atoms1 = np.array(harmonic_atoms1, dtype=np.integer).ravel()
        harmonic_atoms2 = np.array(harmonic_atoms2, dtype=np.integer).ravel()
        
        
        for i, j in izip(harmonic_atoms1, harmonic_atoms2):
            self.draw_bonds.append((i,j))
        
        # do attractive part
        lj_atoms = []
        lj_atoms += parser.get_atom_indices(1, EDGE2_TYPE)
        lj_atoms += parser.get_atom_indices(1, EDGE3_TYPE)
        lj_atoms += parser.get_atom_indices(2, EDGE2_TYPE)
        lj_atoms += parser.get_atom_indices(2, EDGE3_TYPE)
        lj_atoms += parser.get_atom_indices(3, EDGE2_TYPE)
        lj_atoms += parser.get_atom_indices(3, EDGE3_TYPE)
        lj_atoms = np.array(sorted(lj_atoms)).ravel()
        
        self.lj_atoms = lj_atoms
        
        self.extra_atoms = []
        for i in xrange(self.nrigid):
            self.extra_atoms += parser.get_atom_indices(i, OTHER_TYPE)
        
        plate_pot = PlatePotential(harmonic_atoms1, harmonic_atoms2, lj_atoms, k=10)
        # wrap it so it can be used with angle axis coordinates
        pot = RBPotentialWrapper(self.aatopology.cpp_topology, plate_pot)
#            self.aasystem.set_cpp_topology(self.pot.topology)
#            raise Exception
        return pot
    
    def get_potential(self):
        """construct the rigid body potential"""
        try:
            return self.pot
        except AttributeError:
            return self._create_potential()
    
    def get_mindist(self, **kwargs):
        from pele.angleaxis import MinPermDistAACluster
        from pele.angleaxis import TransformAngleAxisCluster, MeasureAngleAxisCluster
        transform = TransformAngleAxisCluster(self.aatopology)
        measure = MeasureAngleAxisCluster(self.aatopology, transform=transform,
                                          permlist=[])
        return MinPermDistAACluster(self.aasystem, measure=measure, transform=transform, 
                                    accuracy=0.1, **kwargs)

    
    def draw(self, rbcoords, index, shift_com=True): # pragma: no cover
        from pele.systems._opengl_tools import draw_atoms, draw_cylinder
        from matplotlib.colors import cnames, hex2color
        
        coords = self.aasystem.to_atomistic(rbcoords).copy()
        coords = coords.reshape([-1,3])
        if shift_com:
            com=np.mean(coords, axis=0)
            coords -= com[np.newaxis, :]
        else:
            com = np.zeros(3)
            
        radius = 0.42

        red = hex2color(cnames["red"])
        c2 = hex2color(cnames["grey"])
        draw_atoms(coords, self.harmonic_atoms, c2, radius=radius)
        draw_atoms(coords, self.lj_atoms, red, radius=radius)
        draw_atoms(coords, self.extra_atoms, c2, radius=radius)

        for i1, i2 in self.draw_bonds:
            draw_cylinder(coords[i1,:], coords[i2,:], color=c2)


def test_bh():
    np.random.seed(0)
    nmol = 4
    system = PlateFolder(nmol)
    db = system.create_database()
    bh = system.get_basinhopping(db)
    bh.run(100)
    m1 = db.minima()[0]
    print m1.coords
    for x in m1.coords:
        print "%.12f," % x,
    print ""
    print m1.energy
    

def test_gui():
    from pele.gui import run_gui
    nmol = 4
    system = PlateFolder(nmol)
    db = system.create_database("tetrahedra.sqlite")
    run_gui(system, db)
    
if __name__ == "__main__":
    test_gui()
#    test_bh()

