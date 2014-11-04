from itertools import izip

import numpy as np
from numpy import cos, sin

from pele.angleaxis import RBTopology, RBSystem, RigidFragment, RBPotentialWrapper
from pele.potentials import BasePotential
from pele.utils import rotations
from plate_potential import PlatePotential

EDGE1_TYPE = "O" 
EDGE2_TYPE = "C"
EDGE3_TYPE = "N"
OTHER_TYPE = "H"

def draw(coords):
    import pele.utils.pymolwrapper as pym
    pym.start()
    pym.draw_spheres(coords, "A", 1)

class HarmonicPotential(BasePotential):
    def __init__(self, atoms1, atoms2):
        self.atoms1 = np.array(atoms1)
        self.atoms2 = np.array(atoms2)
    
    def getEnergy(self, x):
        e, g = self.getEnergyGradient(x)
        return e
    
    def getEnergyGradient(self, x):
        x = x.reshape([-1,3])
        grad = np.zeros(x.shape)
        etot = 0.
        
        for a1, a2 in izip(self.atoms1, self.atoms2):
            dx = x[a1,:] - x[a2,:]
            r2 = np.sum(dx**2)
            etot += 0.5 * r2
            grad[a1,:] += dx 
            grad[a2,:] -= dx
        return etot, grad.ravel()
            
        
class MolAtomIndexParser(object):
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
        return atoms
        
        
        

class CombinePotential(BasePotential):
    def __init__(self, potentials):
        self.potentials = potentials
    
    def getEnergy(self, coords):
        e = 0
        for pot in self.potentials:
            e += pot.getEnergy(coords)
        return e
    
    def getEnergyGradient(self, coords):
        etot = 0
        grad = np.zeros(coords.size)
        for pot in self.potentials:
            e, g = pot.getEnergyGradient(coords)
            etot += e
            grad += g
        return etot, grad.ravel()

def make_plate():
    theta = 60. * np.pi / 180.
    v1 = np.array([1,0,0])
    v2 = np.array([0.5, np.sin(theta), 0])
    
    coords = []
    
    atomi = 0
    edge1 = []
    edge2 = []
    edge3 = []
    emax = 10
    plate = RigidFragment()

    for i in xrange(emax-1):
        for j in xrange(emax-1):
            if i + j >= emax-1: 
                break
            xnew = v1*i + v2*j
            coords.append(xnew)
            if (i == 0 and j == 0 or 
                i == 0 and j == emax-2 or
                i == emax-2 and j == 0):
                atomtype = OTHER_TYPE
            elif i == 0:
                edge1.append(atomi)
                atomtype = EDGE1_TYPE
            elif j == 0:
                edge1.append(atomi)
                atomtype = EDGE2_TYPE
            elif i + j == emax-2:
                edge3.append(atomi)
                atomtype = EDGE3_TYPE
            else:
                atomtype = OTHER_TYPE
            plate.add_atom(atomtype, xnew, 1)
            atomi += 1

    edge1 = np.array(edge1)
    edge2 = np.array(edge2)

    coords = np.array(coords).reshape(-1)
#    draw(coords)
    plate.finalize_setup()
    return plate


class PlateFolder(RBSystem):
    """
    This will build a system class for an OTP (Ortho Ter Phenyl) cluster
    
    OTP is a very simple rigid body molecule defined as 3 Lennard-Jones particles
    connected in a rigid isocolese triangle
    """
    def __init__(self, nmol):
        self.nrigid = nmol
        super(PlateFolder, self).__init__()
        
        self.setup_params(self.params)

    def get_random_configuration(self):
        # js850> this is a bit sketchy because nrigid might not be defined here.
        # probably we can get the number of molecules some other way.
        coords = 100.*np.random.random(6*self.nrigid)
        ca = self.aasystem.coords_adapter(coords)
        for p in ca.rotRigid:
            p[:] = rotations.random_aa()
        return coords


    def setup_aatopology(self):
        """this sets up the topology for the whole rigid body system"""
        topology = RBTopology()
        topology.add_sites([make_plate() for i in xrange(self.nrigid)])
        
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
    
    def get_potential(self):
        """construct the rigid body potential"""
        try:
            return self.pot
        except AttributeError:
            # construct the potential which will compute the energy and gradient in atomistic (cartesian) coordinates
            atomtypes = self.aatopology.get_atomtypes()
            atomtypes = np.array(atomtypes)
            oatoms = np.where(atomtypes == EDGE1_TYPE)[0]
            print atomtypes
            print oatoms
            potentials = []
            parser = MolAtomIndexParser(self.aatopology, self.nrigid)
            
            
            assert self.nrigid == 4
            # do hinges
            atoms1 = []
            atoms2 = []
            e1 = parser.get_atom_indices(0, EDGE1_TYPE)
            e2 = parser.get_atom_indices(1, EDGE1_TYPE)
            atoms1 += list(e1)
            atoms2 += list(e2)
            potentials.append(HarmonicPotential(e1, e2))
             
            e1 = parser.get_atom_indices(0, EDGE2_TYPE)
            e2 = parser.get_atom_indices(2, EDGE1_TYPE)
            atoms1 += list(e1)
            atoms2 += list(e2)
            potentials.append(HarmonicPotential(e1, e2))
            
            e1 = parser.get_atom_indices(0, EDGE3_TYPE)
            e2 = parser.get_atom_indices(3, EDGE1_TYPE)
            atoms1 += list(e1)
            atoms2 += list(e2)
            print "e1", e1
            print "e2", e2
            potentials.append(HarmonicPotential(e1, e2))
            atoms1 = np.array(atoms1, dtype=np.integer).ravel()
            atoms2 = np.array(atoms2, dtype=np.integer).ravel()
            
            # do attractive part
            lj_atoms = []
            lj_atoms += list(parser.get_atom_indices(1, EDGE2_TYPE))
            lj_atoms += list(parser.get_atom_indices(1, EDGE3_TYPE))
            lj_atoms += list(parser.get_atom_indices(2, EDGE2_TYPE))
            lj_atoms += list(parser.get_atom_indices(2, EDGE3_TYPE))
            lj_atoms += list(parser.get_atom_indices(3, EDGE2_TYPE))
            lj_atoms += list(parser.get_atom_indices(3, EDGE3_TYPE))
            print lj_atoms
            lj_atoms = np.array(lj_atoms)
            print lj_atoms
            lj_atoms = lj_atoms.flatten()
            from pele.potentials._lj_cpp import LJCutAtomList
            potentials.append(LJCutAtomList(lj_atoms, rcut=100.))
            
            
#             from pele.potentials._lj_cpp import LJCutAtomList
#             lj = LJCutAtomList(oatoms, rcut=100.)
#             potentials.append(lj)
            from pele.potentials._wca_cpp import WCAAtomList, WCA
#             wca = WCAAtomList(np.array(catoms))
            wca_all = WCA()
            potentials.append(wca_all)
            
#             atomtypes = atomtypes.reshape([self.nrigid,-1])
#             mol1 = atomtypes[0,:]
#             mol2 = atomtypes[1,:]
#             mol1_edge = np.where(mol1 == EDGE2_TYPE)[0]
#             mol2_edge = np.where(mol2 == EDGE2_TYPE)[0] + mol1.size
#             harmonic = HarmonicPotential(mol1_edge, mol2_edge)

            combined_pot = CombinePotential(potentials)
            plate_pot = PlatePotential(atoms1, atoms2, lj_atoms, k=100)
            # wrap it so it can be used with angle axis coordinates
            self.pot = RBPotentialWrapper(self.aatopology.cpp_topology, plate_pot)
#            self.aasystem.set_cpp_topology(self.pot.topology)
#            raise Exception
            return self.pot

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
    
    run_gui(system)
    
if __name__ == "__main__":
    test_gui()
#    test_bh()

