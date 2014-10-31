import numpy as np
from numpy import cos, sin

def draw(coords):
    import pele.utils.pymolwrapper as pym
    pym.start()
    pym.draw_spheres(coords, "A", 1)


from pele.potentials import LJ
from pele.angleaxis import RBTopology, RBSystem, RigidFragment, RBPotentialWrapper


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
            print i, j, i + j
            xnew = v1*i + v2*j
            coords.append(xnew)
            if i == 0:
                edge1.append(atomi)
                atomtype = "O"
            elif j == 0:
                edge1.append(atomi)
                atomtype = "C"
            else:
                atomtype = "N"
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

    def setup_aatopology(self):
        """this sets up the topology for the whole rigid body system"""
        topology = RBTopology()
        topology.add_sites([make_plate() for i in xrange(self.nrigid)])
        
        self.render_scale = 0.2
        self.atom_types = topology.get_atomtypes()
        
        self.draw_bonds = []
        for i in xrange(self.nrigid):
            self.draw_bonds.append((3*i, 3*i+1))
            self.draw_bonds.append((3*i, 3*i+2))
        
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
    
    def get_potential(self):
        """construct the rigid body potential"""
        try:
            return self.pot
        except AttributeError:
            # construct the potential which will compute the energy and gradient in atomistic (cartesian) coordinates
            atomtypes = self.aatopology.get_atomtypes()
            atomtypes = np.array(atomtypes)
            oatoms = np.where(atomtypes == "O")[0]
            catoms = np.where(atomtypes == "C")[0]
            print atomtypes
            print oatoms, catoms
            lj = LJ()
#            from pele.potentials._lj_cpp import LJCutAtomList
            from pele.potentials._wca_cpp import WCAAtomList
            wca = WCAAtomList(np.array(catoms))
            # wrap it so it can be used with angle axis coordinates
            self.pot = RBPotentialWrapper(self.aatopology.cpp_topology, wca)
#            self.aasystem.set_cpp_topology(self.pot.topology)
#            raise Exception
            return self.pot

def test_bh():
    np.random.seed(0)
    nmol = 30
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
    nmol = 2
    system = PlateFolder(nmol)
    
    run_gui(system)
    
if __name__ == "__main__":
    test_gui()
#    test_bh()

