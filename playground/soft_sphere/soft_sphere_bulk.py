import numpy as np

from pele.systems import BaseSystem
from pele.mindist.periodic_exact_match import ExactMatchPeriodic, MeasurePeriodic
from pele.mindist.periodic_mindist import MinDistBulk
from pele.systems.morse_bulk import put_in_box
from pele.potentials import InversePower
from pele.takestep import RandomDisplacement
from pele.transition_states._interpolate import InterpolateLinearMeasure


class SoftSphereSystem(BaseSystem):
    """Binary Lennard Jones potential with periodic boundary conditions"""

    def __init__(self, radii, boxvec, power=2):
        super(SoftSphereSystem, self).__init__(self)

        self.radii = np.array(radii)
        self.natoms = self.radii.size 
        self.boxvec = np.array(boxvec)
        self.periodic = True
        self.potential_kwargs = dict()
        self.power = power
        self.eps = 1.
        
        self.set_params(self.params)
        
        self.params.double_ended_connect.local_connect_params.NEBparams.interpolator = InterpolateLinearMeasure(MeasurePeriodic(self.boxvec))

    def set_params(self, params):
        params.database.accuracy = 1e-4
        params.structural_quench_params.tol = 1e-9
        
#        params.double_ended_connect.local_connect_params.tsSearchParams.iprint = 1
        params.double_ended_connect.local_connect_params.tsSearchParams.hessian_diagonalization = True
        
        params.takestep.stepsize = .4
        
    def get_potential(self):
        return InversePower(self.power, self.eps, self.radii, boxvec=self.boxvec, **self.potential_kwargs)

    def get_random_configuration(self):
        x = np.zeros([self.natoms, self.boxvec.size])
        for i in range(self.boxvec.size):
            x[:, i] = np.random.uniform(-self.boxvec[i] / 2., self.boxvec[i] / 2., self.natoms)
        return x.ravel()

    def draw(self, coordslinear, index):
        from pele.systems._opengl_tools import draw_atomic_binary_polydisperse, draw_box
        put_in_box(coordslinear, self.boxvec)
        draw_atomic_binary_polydisperse(coordslinear, index, bdim=self.boxvec.size, radii=self.radii)
        draw_box(self.boxvec)

    def get_permlist(self):
        return []

    def get_mindist(self):
        measure = MeasurePeriodic(self.boxvec, permlist=self.get_permlist())
        return MinDistBulk(self.boxvec, measure)
    
    def get_compare_exact(self):
        accuracy = 1e-2
        measure = MeasurePeriodic(self.boxvec, self.get_permlist())
        compare = ExactMatchPeriodic(measure, accuracy=accuracy)
        return compare

    def get_orthogonalize_to_zero_eigenvectors(self):
        # TODO: there are some zero eigenvectors which can be removed 
        return None
    
    def get_system_properties(self):
        return dict(natoms=int(self.natoms),
                    radii=self.radii,
                    power=self.power,
                    eps=self.eps,
                    boxvec=self.boxvec,
                    potential="Soft Sphere Bulk",
                    potential_kwargs=self.potential_kwargs,
        )
    
    def get_metric_tensor(self, coords):
        return None

    def get_takestep(self, **kwargs):
        """return the takestep object for use in basinhopping, etc.
        """
        d = dict(self.params.takestep)
        d.update(kwargs)
        kwargs = d
        try:
            stepsize = kwargs.pop("stepsize")
        except KeyError:
            stepsize = 0.6
        takeStep = RandomDisplacement(stepsize=stepsize)
        return takeStep
    

def create_soft_sphere_system_from_db(dbname):
    from pele.storage import Database
    db = Database(dbname, createdb=False)
    
    radii = db.get_property("radii").value()
    boxvec = db.get_property("boxvec").value()
    power = db.get_property("power").value()
    print radii
    
    system = SoftSphereSystem(radii, boxvec, power=power)
    db = system.create_database(dbname, createdb=False)
    
    return system, db


#
# testing only below here
#


def rungui():  # pragma: no cover
    import os
    from pele.gui import run_gui
    dbfname = "test24.sqlite"
    
    if os.path.isfile(dbfname):
        system, db = create_soft_sphere_system_from_db(dbfname) 
    else:
        natoms = 24
        boxl = 3
        boxvec = np.ones(3) * boxl
        # system = MorseCluster(natoms, rho=1.6047, r0=2.8970, A=0.7102, rcut=9.5)
        radii = np.ones(natoms) * .6
        radii += np.random.uniform(-1,1,radii.size) * 1e-1
        system = SoftSphereSystem(radii, boxvec, power=2.5)
        db = system.create_database("test24.sqlite")
    run_gui(system, db)

def plot_potential():
    from matplotlib import pyplot as plt
    natoms = 3
    boxl = 10.
    boxvec = np.ones(3) * boxl
    # system = MorseCluster(natoms, rho=1.6047, r0=2.8970, A=0.7102, rcut=9.5)
    radii = np.ones(natoms) * 1.4
    system = SoftSphereSystem(radii, boxvec, power=4)
    pot = system.get_potential()

    rlist = np.linspace(0,1.5,400)    
    elist = [pot.getEnergy(np.array([0,0.,0.,r,.0,0, 2.5, 0, 0])) for r in rlist]
    plt.plot(rlist, elist)
    print elist
    plt.show()

def test_exact_match():
    natoms = 24
    boxl = 3
    boxvec = np.ones(3) * boxl
    # system = MorseCluster(natoms, rho=1.6047, r0=2.8970, A=0.7102, rcut=9.5)
    radii = np.ones(natoms) * .6
    system = SoftSphereSystem(radii, boxvec, power=2.5)

    x1 = np.genfromtxt("coords1")
    x2 = np.genfromtxt("coords2")

    mindist = system.get_mindist()
    match = system.get_compare_exact()
    
    dist = mindist(x1, x2)[0]
    print dist
    
    ret = match(x1, x2)
    print ret

def test_script():
    natoms = 24
    boxl = 3
    boxvec = np.ones(3) * boxl
    # system = MorseCluster(natoms, rho=1.6047, r0=2.8970, A=0.7102, rcut=9.5)
    radii = np.ones(natoms) * .6
    system = SoftSphereSystem(radii, boxvec, power=2.5)
    
    db = system.create_database()
    
    bh = system.get_basinhopping(db)
    bh.run(100)
    
    from pele.landscape import ConnectManager
    manager = ConnectManager(db, strategy="random")
    for i in xrange(10):
        try:
            m1, m2 = manager.get_connect_job()
        except manager.NoMoreConnectionsError:
            break
        
        connect = system.get_double_ended_connect(m1, m2, db)
        connect.connect()
            
def print_info(dbfname="ss60.sqlite"):
    system, db = create_soft_sphere_system_from_db(dbfname)
    
    import networkx as nx
    from pele.landscape import database2graph
    graph = database2graph(db)
    
    cc = nx.connected_components(graph)
    c = cc[0]
    for i in xrange(3):
        print c[i].energy
      

if __name__ == "__main__":
#    plot_potential()
#    rungui()
#    test_exact_match()
#    test_script()
    print_info()
