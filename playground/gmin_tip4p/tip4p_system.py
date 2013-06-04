import numpy as np
from pele.potentials import GMINPotential
import gmin_ as GMIN
from copy import deepcopy
import tip4p

from pele.angleaxis import RBTopology, RBSystem

class TIP4PSystem(RBSystem):
    def __init__(self):
        RBSystem.__init__(self)
    
    def setup_aatopology(self):
        GMIN.initialize()
        pot = GMINPotential(GMIN)
        coords = pot.getCoords()        
        nrigid = coords.size / 6

        print "I have %d water molecules in the system"%nrigid
        print "The initial energy is", pot.getEnergy(coords)

        water = tip4p.water()
        
        system = RBTopology()
        system.add_sites([deepcopy(water) for i in xrange(nrigid)])
        self.potential = pot
        self.nrigid = nrigid
        
        self.render_scale = 0.3
        self.atom_types = system.get_atomtypes()
        
        self.draw_bonds = []
        for i in xrange(nrigid):
            self.draw_bonds.append((3*i, 3*i+1))
            self.draw_bonds.append((3*i, 3*i+2))
    
        return system
    
    def get_potential(self):
        return self.potential
    
if __name__ == "__main__":
    import pele.gui.run as gr
    gr.run_gui(TIP4PSystem, db="tip4p_8.sqlite")
