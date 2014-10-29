import numpy as np
from numpy import cos, sin, pi

#import gmin_ as GMIN

#from pele.potentials import LJ
from pele.angleaxis import RBTopologyBulk, RBSystem, RigidFragmentBulk, RBPotentialWrapper
from pele.potentials.ljcut import LJCut
from pele.mindist.minpermdist_stochastic import MinPermDistBulk
from pele.mindist.periodic_exact_match import MeasurePeriodic

def put_in_box(x, boxvec):
    x = x.reshape(-1, boxvec.size)
    x -= boxvec * np.round(x / boxvec)

class OTPBulk(RBSystem):
    """
    This will build a system class for a system of OTP (Ortho Ter Phenyl) molecules
    
    OTP is a very simple rigid body molecule described as 3 Lennard-Jones particles
    connected in a rigid isosceles triangle
    """
    def __init__(self, nmol,boxl,rcut):
        self.nrigid = nmol
        self.boxvec=boxl
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
        topology = RBTopologyBulk(self.boxvec)   # sn402: changed
        #topology = RBTopology()       
        topology.add_sites([self.make_otp() for i in xrange(self.nrigid)])
        
        self.render_scale = 0.2
        self.atom_types = topology.get_atomtypes()
        
        self.draw_bonds = []
        for i in xrange(self.nrigid):
            self.draw_bonds.append((3*i, 3*i+1))
            self.draw_bonds.append((3*i, 3*i+2))
        
        topology.finalize_setup()
    
        return topology
    
    def get_random_configuration(self):
        x = np.zeros([self.nrigid,6])
        for i in range(3):
            x[:,i] = np.random.uniform(-self.boxvec[i]/2., self.boxvec[i]/2., self.nrigid)
        for i in range(3,6):
            x[:,i] = 5.*np.random.random(self.nrigid)
        return x.flatten()

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
            # construct the potential which will compute the energy and gradient in atomistic (cartesian) coordinates
            cartesian_potential = LJCut(rcut=self.cut, boxl=self.boxvec)
            # wrap it so it can be used with angle axis coordinates
            self.pot = RBPotentialWrapper(self.aatopology.cpp_topology, cartesian_potential)
#            self.aasystem.set_cpp_topology(self.pot.topology)
            return self.pot
        
    def get_mindist(self, **kwargs):
        measure = MeasurePeriodic(self.boxvec)
        return MinPermDistBulk(self.boxvec, measure, niter=10, verbose=False, tol=0.01, 
                 accuracy=0.01)

        

def test_bh():
    np.random.seed(0)
    nmol = 5
    boxl = np.array([15,10,5])
    rcut = 2.5
    system = OTPBulk(nmol,boxl,rcut)   
    db = system.create_database()
    bh = system.get_basinhopping(db)   # sn402: just overload get_random_configuration() in 
    # this file when it's time to specify the initial coordinates of OTP.
    bh.run(100)
    m1 = db.minima()[0]
    print m1.coords
    for x in m1.coords:
        print "%.12f," % x,
    print ""
    print m1.energy
    print db.minima()[1].energy
    print db.minima()[2].energy   
    return db

def test_gui():
    from pele.gui import run_gui
    nmol = 5
    system = OTPBulk(nmol,np.array([15,10,5]),2.5)
    
    run_gui(system)
    
def test_PBCs():
    np.random.seed(0)
    nmol = 2
    boxl = np.array([5,5,5])
    rcut = 2.5
    system = OTPBulk(nmol,boxl,rcut)
    #coords1 = system.get_random_configuration()
    #coords2 = system.get_random_configuration() 
    coords1 = np.array([-2.,1.,1.,1.,1.,1.,0.,0.,0.,0.,0.,0.]) 
    coords2 = np.array([3,1.,1.,1.,1.,1.,0.,0.,0.,0.,0.,0.])  
    print coords1
    print coords2
    a = system.aatopology.distance_squared(coords1,coords2) 
    # Note, for displacements without rotation this gives 3x the square displacement of the com: 
    # there are 3 atoms being displaced.
    print a
    
def test_mindist():
    nmol = 2
    boxl = np.array([5,5,5])
    rcut = 2.5
    system = OTPBulk(nmol,boxl,rcut) 
    #print system.aatopology.boxvec  
    coords1 = np.array([0.,0.,0.,0.,-2.,-2.,0.,0.,0.,0.,0.,0.])
    coords2 = np.array([1.,0.,0.,0.,2.,2.,0.,0.,0.,0.,0.,0.])    
    #print coords
    
    import pele.mindist.periodic_exact_match as pd
    a = pd.MeasurePeriodicRigid(boxl, system.aatopology)
    b = a.get_dist(coords1, coords2)
    print b
    
def test_connect():

    nmol = 5
    boxl = np.array([10,10,5])
    rcut = 2.5
    system = OTPBulk(nmol,boxl,rcut)   

    db = test_bh()
    X1 = db.minima()[0].coords
    X2 = db.minima()[1].coords
    
    import pele.mindist.periodic_exact_match as md
    a = md.MeasurePeriodicRigid(boxl,system.aatopology)
    b = MinPermDistBulk(boxl, a)
      
    dist, x1, x2 = b(X1,X2)
    
    min1, min2 = db.minima()[0], db.minima()[1]
    #from pele.landscape import ConnectManager
    #manager = ConnectManager(db, strategy="gmin")
    #for i in xrange(db.number_of_minima()-1):
    #    min1, min2 = manager.get_connect_job()
    connect = system.get_double_ended_connect(min1, min2, db)
    connect.connect()
        
    #from pele.utils.disconnectivity_graph import DisconnectivityGraph, database2graph
    #import matplotlib.pyplot as plt
    #convert the database to a networkx graph
    #graph = database2graph(db)
    #dg = DisconnectivityGraph(graph, nlevels=3, center_gmin=True)
    #dg.calculate()
    #dg.plot()
    #plt.show()
     
if __name__ == "__main__":
#    test_gui()
#    test_bh()
    test_connect()
#    test_PBCs()
#    test_mindist()
