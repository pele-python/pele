import tempfile
from copy import deepcopy

import numpy as np
from numpy import cos, sin, pi

import gmin_ as GMIN

from pygmin.potentials import GMINPotential
from pygmin.potentials import LJ
from pygmin.angleaxis import RBTopology, RBSystem, RigidFragment
from pygmin.utils.xyz import write_xyz


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
#        super(OTPSystem, self)
        self.nmol = nmol
        if boxl is None:
            self.periodic = False
            self.boxl = 1e100
        else:
            self.periodic = True
            self.boxl = boxl
        super(OTPSystem, self).__init__()
        
        self.setup_params(self.params)

    def get_rigid_fragment(self):
        return RigidFragment()

    def make_otp(self):
        otp = self.get_rigid_fragment()
        otp.add_atom("O", np.array([0.0, -2./3 * np.sin( 7.*pi/24.), 0.0]), 1.)
        otp.add_atom("O", np.array([cos( 7.*pi/24.),  1./3. * sin( 7.* pi/24.), 0.0]), 1.)
        otp.add_atom("O", np.array([-cos( 7.* pi/24.),  1./3. * sin( 7.*pi/24), 0.0]), 1.)
        otp.finalize_setup()
        return otp

    
    def setup_params(self, params):
        params.gui["basinhopping_nsteps"] = 1000
        
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
        otp = self.make_otp()
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

    def load_coords_pymol(self, coordslist, oname, index=1):
        """load the coords into pymol
        
        the new object must be named oname so we can manipulate it later
                        
        Parameters
        ----------
        coordslist : list of arrays
        oname : str
            the new pymol object must be named oname so it can be manipulated
            later
        index : int
            we can have more than one molecule on the screen at one time.  index tells
            which one to draw.  They are viewed at the same time, so should be
            visually distinct, e.g. different colors.  accepted values are 1 or 2
        
        Notes
        -----
        the implementation here is a bit hacky.  we create a temporary xyz file from coords
        and load the molecule in pymol from this file.  
        """
        # pymol is imported here so you can do, e.g. basinhopping without installing pymol
        import pymol 

        # create the temporary file
        suffix = ".xyz"
        f = tempfile.NamedTemporaryFile(mode="w", suffix=suffix)
        fname = f.name
                
        # write the coords into the xyz file
        from pygmin.mindist import CoMToOrigin
        for coords in coordslist:
            if hasattr(self, "atom_types"):
                atom_types = self.atom_types
            else:
                atom_types = ["O"]
            atom_coords = self.aasystem.to_atomistic(coords)
#            atom_coords = CoMToOrigin(coords.copy())
            write_xyz(f, atom_coords, title=oname, atomtypes=atom_types)#["C", "N", "N"])
        f.flush()
                
        # load the molecule from the temporary file
        pymol.cmd.load(fname)
        
        # get name of the object just create and change it to oname
        objects = pymol.cmd.get_object_list()
        objectname = objects[-1]
        pymol.cmd.set_name(objectname, oname)
        
        #set the representation
        pymol.cmd.hide("everything", oname)
        pymol.cmd.show("spheres", oname)
        
        
        pymol.cmd.unbond(oname, oname)
        
        if hasattr(self, "draw_bonds"):
            for i1, i2 in self.draw_bonds:
                pymol.cmd.bond("id "+str(i1+1), "id "+str(i2+1))
            pymol.cmd.show("lines", oname)

        # set the color according to index
        if index == 1:
            pymol.cmd.color("red", oname)
        else:
            pymol.cmd.color("gray", oname)

        pymol.cmd.set("sphere_scale", value=0.2, selection=oname)

def rungui(system, db=None):
    import pygmin.gui.run as gr
    from pygmin.storage import Database
    gr.run_gui(system, db=db)
    
if __name__ == "__main__":
    nmol = 7
    periodic = False
    if periodic:
        boxl = 10.
        dbname = "otp_N%d_boxl%.2f.db"
    else:
        dbname = "otp_N%d.db" % nmol
        boxl = None
    system = OTPSystem(nmol, boxl=boxl)
#    bh = system.get_basinhopping()
#    bh.run(10)
    
    rungui(system, db=dbname)
    
