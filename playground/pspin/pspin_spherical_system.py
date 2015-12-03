from __future__ import division
import numpy as np

from pele.potentials import MeanFieldPSpinSpherical
from pele.systems import BaseSystem
from pele.landscape import smooth_path
from scipy.misc import factorial

def normalize_spins(x):
    x /= (np.linalg.norm(x)/np.sqrt(len(x)))
    return x

def spin_distance_1d(x1, x2):
    x1 = normalize_spins(x1)
    x2 = normalize_spins(x2)
    return np.linalg.norm(x1 - x2)

def spin_mindist_1d(x1, x2):
    x1 = normalize_spins(x1)
    x2 = normalize_spins(x2)
    return np.linalg.norm(x1-x2), x1, x2

class MeanFieldPSpinSphericalSystem(BaseSystem):
    def __init__(self, nspins, p=3):
        BaseSystem.__init__(self)
        self.nspins = nspins
        self.p = p
        self.interactions = self.get_interactions(self.nspins, self.p)
        self.pot = self.get_potential()
               
        self.setup_params(self.params)

    def setup_params(self, params):
#        params.takestep.stepsize = np.pi# / 2.
        params.takestep.verbose = True
#        params.double_ended_connect.local_connect_params.NEBparams.interpolator = interpolate_spins
        params.double_ended_connect.local_connect_params.NEBparams.image_density = .8
        params.double_ended_connect.local_connect_params.NEBparams.iter_density = 50.
        params.double_ended_connect.local_connect_params.NEBparams.reinterpolate = 50
        params.double_ended_connect.local_connect_params.NEBparams.adaptive_nimages = True
        params.double_ended_connect.local_connect_params.NEBparams.adaptive_niter = False
#        params.double_ended_connect.local_connect_params.NEBparams.distance = spin3d_distance
        params.structural_quench_params.tol = 1e-6
        params.database.overwrite_properties = False
        
        params.basinhopping.insert_rejected = True

    def get_system_properties(self):
        return dict(potential="PSpinSPherical model",
                    nspins=self.nspins,
                    p=self.p,
                    interactions=self.interactions,
                    )
        
    def get_interactions(self, nspins, p):
        assert p == 3
        interactions = np.random.normal(0, np.sqrt(factorial(p)), [nspins for i in xrange(p)])
        for i in xrange(nspins):
            for j in xrange(i, nspins):
                for k in xrange(j, nspins):
                    interactions[k][i][j] = interactions[i][j][k]
                    interactions[k][j][i] = interactions[i][j][k]
                    interactions[j][k][i] = interactions[i][j][k]
                    interactions[i][k][j] = interactions[i][j][k]
                    interactions[j][i][k] = interactions[i][j][k]
        return interactions.flatten()

    def get_potential(self, tol=1e-6):
        try:
            return self.pot
        except AttributeError:
            self.pot = MeanFieldPSpinSpherical(self.interactions, self.nspins, self.p, tol=tol)
            return self.pot

    def get_orthogonalize_to_zero_eigenvectors(self):
        return None
    
    def get_metric_tensor(self, coords):
        return None
    
    def get_nzero_modes(self):
        return 0

    def get_pgorder(self, coords):
        return 1
    
    def get_mindist(self):
        return spin_mindist_1d

    def get_compare_exact(self):
        """
        are they the same minima?
        """
        mindist = self.get_mindist()
        return lambda x1, x2: mindist(x1, x2)[0]/np.sqrt(self.nspins) < 1e-3

    def smooth_path(self, path, **kwargs):
        mindist = self.get_mindist()
        return smooth_path(path, mindist, **kwargs)

    def get_random_configuration(self):
        coords = np.random.normal(0, 1, self.nspins)
        coords /= (np.linalg.norm(coords)/np.sqrt(self.nspins))
        return coords
    
    def create_database(self, *args, **kwargs):
        return BaseSystem.create_database(self, *args, **kwargs)

    def get_takestep(self, **kwargs):
        """return the takestep object for use in basinhopping, etc.
        
        default is random displacement with adaptive step size 
        adaptive temperature
        
        See Also
        --------
        pele.takestep
        """
        #return super(MeanFieldPSpinSphericalSystem, self).get_takestep(**kwargs)
        # if no disorder, turn off adaptive step and temperature.
        from pele.takestep import RandomDisplacement
        kwargs = dict(self.params["takestep"].items() + kwargs.items())
        try:
            stepsize = kwargs.pop("stepsize")
        except KeyError:
            stepsize = self.nspins
        takeStep = RandomDisplacement(stepsize=stepsize)
        return takeStep

    def draw(self, coords, index):
        pass
    
def normalize_spins_db(db):
    for m in db.minima():
        x = normalize_spins(m.coords)
        print np.max(x), np.min(x)
        m.coords = x
    db.session.commit()
    
def run_gui():
    from pele.gui import run_gui
    system = MeanFieldPSpinSphericalSystem(20, p=3)
    run_gui(system)

#def run_gui_db(dbname="xy_10x10.sqlite"):
#    from pele.gui import run_gui
#    from pele.storage import Database
#    try:
#        db = Database(dbname, createdb=False)
#        phases = db.get_property("phases").value()
#    except IOError:
#        phases=None
#    system = XYModlelSystem(dim=[10,10], phi_disorder=np.pi, phases=phases)
#    run_gui(system, db=dbname)

if __name__ == "__main__":
    run_gui()
#    test_potential()
#    from pele.storage import Database
#    db = Database("20x20_no_disorder.sqlite")
#    normalize_spins_db(db)
#    run_gui_db()
#    run_gui_nodisorder()
        
    
    
