import numpy as np
from pygmin.potentials import GMINPotential
import gmin_ as GMIN
from copy import deepcopy

from pygmin.angleaxis import RigidFragment, RBTopology
from pygmin.angleaxis import RBTopology, RBSystem

from math import sin, cos, pi
from pygmin.optimize import fire

def create_pap():
    pap = RigidFragment()
    rho   = 0.3572
    theta = 104.52/180.0*pi      
    pap.add_atom("O", np.array([0., 0., 0.]), 1.)
    pap.add_atom("H", np.array([-0.038490, 0.1204928, 0.3794728]), 1.)
    pap.add_atom("C", np.array([-0.038490, -0.1204928, -0.3794728]), 1.)
    pap.finalize_setup(shift_com=False)
    
    print "inversion:\n", pap.inversion
    print "symmetry:\n", pap.symmetries
    pap.inversion=None
    
    return pap

class PAPSystem(RBSystem):
    def __init__(self):
        RBSystem.__init__(self)
        
        NEBparams = self.params.double_ended_connect.local_connect_params.NEBparams
        self.params.double_ended_connect.local_connect_params.NEBparams.distance=self.aasystem.neb_distance

        NEBparams.max_images=200
        NEBparams.image_density=10.
        NEBparams.iter_density=6
        NEBparams.quenchRoutine = fire
        del NEBparams.NEBquenchParams["maxErise"]
        NEBparams.NEBquenchParams["maxstep"]=0.01
        
        fts = self.params.double_ended_connect.local_connect_params.tsSearchParams
        fts["tol"]=0.01
    
    def setup_aatopology(self):
        GMIN.initialize()
        pot = GMINPotential(GMIN)
        coords = pot.getCoords()        
        nrigid = coords.size / 6

        print "I have %d PAP molecules in the system"%nrigid
        print "The initial energy is", pot.getEnergy(coords)

        water = create_pap()
        
        system = RBTopology()
        system.add_sites([deepcopy(water) for i in xrange(nrigid)])
        self.potential = pot
        self.nrigid = nrigid
        
        self.render_scale = 0.1
        self.atom_types = system.get_atomtypes()
        
        self.draw_bonds = []
        for i in xrange(nrigid):
            self.draw_bonds.append((3*i, 3*i+1))
            self.draw_bonds.append((3*i, 3*i+2))
    
        return system
    
    def get_potential(self):
        return self.potential
    
if __name__ == "__main__":
    import pygmin.gui.run as gr
    gr.run_gui(PAPSystem, db="pap.sqlite")
