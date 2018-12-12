from __future__ import print_function
import numpy as np
from numpy import cos, sin, pi

#import gmin_ as GMIN

#from pele.potentials import LJ
from pele.angleaxis import RBTopologyBulk, RBSystem, RigidFragmentBulk, RBPotentialWrapper
#from pele.potentials.ljcut import LJCut
from pele.potentials._lj_cpp import LJCutCellLists, LJCut
from pele.angleaxis.bulk_rigid_mindist import MinDistBulkRigid
from pele.angleaxis.aaperiodicttransforms import MeasurePeriodicRigid,\
    TransformPeriodicRigid

class OTPBulk(RBSystem):
    """
    This will build a system class for a system of OTP (Ortho Ter Phenyl) molecules
    
    OTP is a very simple rigid body molecule described as 3 Lennard-Jones particles
    connected in a rigid isosceles triangle
    """
    def __init__(self, nmol,boxvec,rcut):
        self.nrigid = nmol
        self.boxvec = np.array(boxvec, dtype=float)
        self.cut=rcut
        
        super(OTPBulk, self).__init__()      
        self.setup_params(self.params)


    def make_otp(self):
        """this constructs a single OTP molecule"""
        otp = RigidFragmentBulk(self.boxvec)   # sn402: changed
        otp.add_atom("O", np.array([0.0, -2./3 * np.sin( 7.*pi/24.), 0.0]), 1.)
        otp.add_atom("O", np.array([cos( 7.*pi/24.),  1./3. * sin( 7.* pi/24.), 0.0]), 1.)
        otp.add_atom("O", np.array([-cos( 7.* pi/24.),  1./3. * sin( 7.*pi/24), 0.0]), 1.)
        otp.finalize_setup()
        return otp
        
    def setup_aatopology(self):
        """this sets up the topology for the whole rigid body system"""
        topology = RBTopologyBulk(self.boxvec)
 
        topology.add_sites([self.make_otp() for i in range(self.nrigid)])
        
        self.render_scale = 0.2
        self.atom_types = topology.get_atomtypes()
        
        self.draw_bonds = []
        for i in range(self.nrigid):
            self.draw_bonds.append((3*i, 3*i+1))
            self.draw_bonds.append((3*i, 3*i+2))
        
        topology.finalize_setup()  
        return topology
    
    def get_random_configuration(self):
        """ Returns an array containing random periodic com/aa coordinates."""
        x = np.zeros([self.nrigid*2,3])
        for i in range(self.nrigid):
            for j in range(3):
                x[i][j] = np.random.uniform(-self.boxvec[j]/2., self.boxvec[j]/2.)
        for i in range(self.nrigid,2*self.nrigid):
            x[i] = 5.*np.random.random(3)
        return x.ravel()
        
    def configuration_from_file(self, fileobj, angleaxis=True):
        """ Returns an array of com/aa coordinates as read in from a file.
        
            Parameters
            ----------
            fileobj - file
            A file object corresponding to the input file
            angleaxis - boolean, optional 
            If this is true, the input file should consist of 3*nrigid 
            centre-of-mass coordinates followed by 3*nrigid angle-axis vector components. 
            The exact shape does not matter.
            Otherwise, the input file consists of 9*nrigid atomistic cartesian coordinates.
        """
        x = []

        for line in fileobj:
            y = line.split()
            for item in y:
                x.append(float(item))

        if angleaxis is True:
            if (len(x)%(6*self.nrigid) != 0):
                raise IOError("Input file is the wrong length: read in {} values", len(x))
            return np.array(x).ravel()
        else:
            if(len(x)%(9*self.nrigid) != 0):
                raise IOError("Input file is the wrong length")
            return np.array

    def setup_params(self, params):
        """set some system dependent parameters to improve algorithm performance"""
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
            # construct the potential which will compute the energy and gradient 
            # in atomistic (cartesian) coordinates
            # NOTE: Currently the LJCut potential only deals with cubic boxes
            cartesian_potential = LJCut(rcut=self.cut, boxvec=self.boxvec)
#            cartesian_potential = LJCutCellLists(rcut=self.cut, boxvec=self.boxvec)
            # wrap it so it can be used with angle axis coordinates
            self.pot = RBPotentialWrapper(self.aatopology.cpp_topology, cartesian_potential)
#            self.aasystem.set_cpp_topology(self.pot.topology)
            return self.pot
        
    def get_mindist(self, **kwargs):
        measure = MeasurePeriodicRigid(self.aatopology, transform=TransformPeriodicRigid())
        return MinDistBulkRigid(self.boxvec, measure, niter=10, transform=TransformPeriodicRigid(), 
                               verbose=False, tol=0.01, accuracy=0.01)

    def draw(self, rbcoords, index):
        from pele.systems._opengl_tools import draw_box
        from pele.systems.morse_bulk import put_in_box
        cc = self.aatopology.coords_adapter(rbcoords)
        put_in_box(cc.posRigid, self.boxvec)
        super(OTPBulk, self).draw(rbcoords, index, shift_com=False)
        draw_box(self.boxvec)

#
# only testing below here
#
        
def test_bh():  # pragma: no cover
    np.random.seed(0)
    nmol = 5
    boxvec = np.array([5.,5.,5.])
    rcut = 2.5
    system = OTPBulk(nmol,boxvec,rcut)   
    db = system.create_database()
    bh = system.get_basinhopping(db)   # sn402: just overload get_random_configuration() in 
    # this file when it's time to specify the initial coordinates of OTP.
    bh.run(100)
    m1 = db.minima()[0]
    print(m1.coords)
    for x in m1.coords:
        print("%.12f," % x, end=' ')
    print("")
    m2 = db.minima()[1]
    print(m2.coords)
    for x in m2.coords:
        print("%.12f," % x, end=' ')
    print("")   
    
    print(m1.energy)
    print(db.minima()[1].energy)
    print(db.minima()[2].energy)      
 

def test_gui():  # pragma: no cover
    from pele.gui import run_gui
    nmol = 20
    boxvec = np.array([6,6,6])
    rcut = 2.5
    system = OTPBulk(nmol, boxvec, rcut)
    
    run_gui(system)
       
def test_mindist():  # pragma: no cover
    nmol = 2
    boxvec = np.array([15,10,5])
    rcut = 2.5
    system = OTPBulk(nmol,boxvec,rcut) 
    #print system.aatopology.boxvec  
    coords1 = np.array([0.,0.,0.,0.,-2.,-2.,0.,0.,0.,0.,0.,0.])
    coords2 = np.array([1.,0.,0.,0.,2.,2.,0.,0.,0.,0.,0.,0.])    
    #print coords
    
    import pele.angleaxis.aaperiodicttransforms as pd
    a = pd.MeasurePeriodicRigid(boxvec, system.aatopology)
    b = a.get_dist(coords1, coords2)
    print(b)
    
def test_connect():  # pragma: no cover

    nmol = 5
    boxvec = np.array([5,5,5])
    rcut = 2.5
    system = OTPBulk(nmol,boxvec,rcut)   

    db = test_bh()
#     X1 = db.minima()[0].coords
#     X2 = db.minima()[1].coords
    
#     import pele.angleaxis.aaperiodicttransforms as md
#     a = md.MeasurePeriodicRigid(system.aatopology, transform=TransformPeriodicRigid())
#     b = MinPermDistBulk(boxvec, a, transform=TransformPeriodicRigid())     
#     dist, x1, x2 = b(X1,X2)
    
    min1, min2 = db.minima()[0], db.minima()[1]
#     from pele.landscape import ConnectManager
#     manager = ConnectManager(db, strategy="gmin")
#     for i in xrange(db.number_of_minima()-1):
#         min1, min2 = manager.get_connect_job()
    connect = system.get_double_ended_connect(min1, min2, db)
    connect.connect()
    from pele.utils.disconnectivity_graph import DisconnectivityGraph, database2graph
    import matplotlib.pyplot as plt
#     convert the database to a networkx graph
    graph = database2graph(db)
    dg = DisconnectivityGraph(graph, nlevels=3, center_gmin=True)
    dg.calculate()
    dg.plot()
    plt.show()          

if __name__ == "__main__":
     test_gui()
#    test_bh()
#     test_connect()
#    test_PBCs()
#    test_mindist()

