import numpy as np
from pygmin import gui
import oxdnagmin_ as GMIN
import pygmin.gui.run as gr
from pygmin.utils.rbtools import CoordsAdapter
from pygmin import takestep
from pygmin.utils import rotations
from pygmin.potentials import GMINPotential
from pygmin.transition_states import NEB, InterpolatedPath
from pygmin.angleaxis import RBSystem
from pygmin.angleaxis.molecules import create_water
from pygmin.angleaxis import RBTopology, RigidFragment
from copy import deepcopy

def create_base():
    base = RigidFragment()
    base.add_atom("O", np.array([-0.4, 0., 0.]), 1.)
    base.add_atom("H", np.array([0.4, 0., 0.]), 1.)
    base.finalize_setup(shift_com=False)
    
    print "inversion:\n", base.inversion
    print "symmetry:\n", base.symmetries
    base.inversion=None
    
    return base


class OXDNASystem(RBSystem):
    
    def setup_aatopology(self):
        GMIN.initialize()
        pot = GMINPotential(GMIN)
        coords = pot.getCoords()        
        nrigid = coords.size / 6

        print "I have %d water molecules in the system"%nrigid
        print "The initial energy is", pot.getEnergy(coords)

        water = create_base()
        
        system = RBTopology()
        system.add_sites([deepcopy(water) for i in xrange(nrigid)])
        self.potential = pot
        self.nrigid = nrigid
        
        self.render_scale = 0.15
        self.atom_types = system.get_atomtypes()
        
        self.draw_bonds = []
        for i in xrange(nrigid-1):
            self.draw_bonds.append((2*i, 2*i+1))
            self.draw_bonds.append((2*i, 2*i+2))
    
        return system

    def get_potential(self):
        return self.potential
        
#    def get_mindist(self):
#        mindist = RBSystem.get_mindist(self)
#        mindist.transform
        

if __name__ == "__main__":
    
    gr.run_gui(OXDNASystem, db="oxdna.sqlite")