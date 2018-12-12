from __future__ import print_function
import numpy as np
from numpy import cos, sin, pi

#import gmin_ as GMIN

from pele.potentials import LJ
from pele.angleaxis import RBTopology, RBSystem, RigidFragment, RBPotentialWrapper

def make_otp():
    """this constructs a single OTP molecule"""
    otp = RigidFragment()
    otp.add_atom("O", np.array([0.0, -2./3 * np.sin( 7.*pi/24.), 0.0]), 1.)
    otp.add_atom("O", np.array([cos( 7.*pi/24.),  1./3. * sin( 7.* pi/24.), 0.0]), 1.)
    otp.add_atom("O", np.array([-cos( 7.* pi/24.),  1./3. * sin( 7.*pi/24), 0.0]), 1.)
    otp.finalize_setup()
    return otp

class OTPCluster(RBSystem):
    """
    This will build a system class for an OTP (Ortho Ter Phenyl) cluster
    
    OTP is a very simple rigid body molecule defined as 3 Lennard-Jones particles
    connected in a rigid isocolese triangle
    """
    def __init__(self, nmol):
        self.nrigid = nmol
        super(OTPCluster, self).__init__()
        
        self.setup_params(self.params)
        
    def setup_aatopology(self):
        """this sets up the topology for the whole rigid body system"""
        topology = RBTopology()
        topology.add_sites([make_otp() for _ in range(self.nrigid)])
        
        self.render_scale = 0.2
        self.atom_types = topology.get_atomtypes()
        
        self.draw_bonds = []
        for i in range(self.nrigid):
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
            cartesian_potential = LJ()
            # wrap it so it can be used with angle axis coordinates
            self.pot = RBPotentialWrapper(self.aatopology.cpp_topology, cartesian_potential)
#            self.aasystem.set_cpp_topology(self.pot.topology)
            return self.pot
    
    def load_coords_pymol(self, *args, **kwargs):
        from . import pymol
        RBSystem.load_coords_pymol(self, *args, **kwargs)
        # draw the spheres slightly smaller
        pymol.cmd.set("sphere_scale", value=.25)        
    

def test_bh():
    np.random.seed(0)
    nmol = 10
    system = OTPCluster(nmol)
    db = system.create_database()
    bh = system.get_basinhopping(db)
    bh.run(100)
    m1 = db.minima()[0]
    print(m1.coords)
    for x in m1.coords:
        print("%.12f," % x, end=' ')
    print("")
    print(m1.energy)
    

def test_gui():
    from pele.gui import run_gui
    nmol = 5
    system = OTPCluster(nmol)
    
    run_gui(system)
    
if __name__ == "__main__":
    test_gui()
#    test_bh()

