import numpy as np
from numpy import cos, sin, pi
from pygmin.potentials import GMINPotential
import gmin_ as GMIN
from copy import deepcopy
from pygmin.potentials import LJ
from pygmin.angleaxis import RBTopology, RBSystem, RigidFragment

def make_otp():
    otp = RigidFragment()
    otp.add_atom("O", np.array([0.0, -2./3 * np.sin( 7.*pi/24.), 0.0]), 1.)
    otp.add_atom("O", np.array([cos( 7.*pi/24.),  1./3. * sin( 7.* pi/24.), 0.0]), 1.)
    otp.add_atom("O", np.array([-cos( 7.* pi/24.),  1./3. * sin( 7.*pi/24), 0.0]), 1.)
    otp.finalize_setup()
    return otp

class OTPPotential(LJ):
    def __init__(self, aatopology):
        self.aatopology = aatopology
        LJ.__init__(self)
        
    def getEnergy(self, coords):
        atom_coords = self.aatopology.to_atomistic(coords)
        return LJ.getEnergy(self, atom_coords.flatten())

    def getEnergyGradient(self, coords):
        atom_coords = self.aatopology.to_atomistic(coords)
        e, atom_grad = LJ.getEnergyGradient(self, atom_coords.flatten())
        grad = self.aatopology.transform_gradient(coords, atom_grad)
        return e, grad

class OTPSystem(RBSystem):
    def __init__(self, nmol, boxl=None):
        self.nmol = nmol
        if boxl is None:
            self.periodic = False
            self.boxl = 1e100
        else:
            self.periodic = True
            self.boxl = boxl
        RBSystem.__init__(self)
    
    def get_random_coordinates(self):
        coords = np.zeros([self.nmol*2, 3])
        coords[:self.nmol,:] = np.random.uniform(-1,1,[self.nmol,3]) * (self.nmol*3)**(-1./3) * 1.5
        from pygmin.utils.rotations import random_aa
        for i in range(self.nmol, self.nmol*2):
            coords[i,:] = random_aa()
        return coords.flatten()
    
    def write_coords_data(self):
        coords = self.get_random_coordinates()
        coords = coords.reshape(-1,3)
        with open("coords", "w") as fout:
            for xyz in coords:
                fout.write( "%f %f %f\n" % tuple(xyz) )
        
        data = "LWOTP 2.614\n"
        if self.periodic:
            data += "periodic %f %f %f\n" % (self.boxl, self.boxl, self.boxl)
        
        with open("data", "w") as fout:
            fout.write(data)
    
    def setup_aatopology(self):
        self.write_coords_data()
        GMIN.initialize()        
        self.pot = GMINPotential(GMIN)
        coords = self.pot.getCoords()
        self.nrigid = coords.size/6
        assert(self.nrigid == self.nmol)
        #self.nrigid = self.nmol
        otp = make_otp()
        topology = RBTopology()
        topology.add_sites([deepcopy(otp) for i in xrange(self.nrigid)])
        
        self.render_scale = 0.2
        self.atom_types = topology.get_atomtypes()
        
        self.draw_bonds = []
        for i in xrange(self.nrigid):
            self.draw_bonds.append((3*i, 3*i+1))
            self.draw_bonds.append((3*i, 3*i+2))
    
        self.params.double_ended_connect.local_connect_params.tsSearchParams.iprint = 10
        return topology
    
    def get_potential(self):
        #return OTPPotential(self.aasystem)
        return self.pot
    
    def __call__(self):
        return self

def rungui(system, db=None):
    import pygmin.gui.run as gr
    gr.run_gui(system, db=db)
    
if __name__ == "__main__":
    nmol = 20
    system = OTPSystem(nmol, boxl=10)
#    bh = system.get_basinhopping()
#    bh.run(10)
    
    rungui(system)
    
