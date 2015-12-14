from __future__ import division
import numpy as np
import cmath
from numba import jit

from pele.potentials import MeanFieldPSpinSpherical
from pele.systems import BaseSystem
from pele.landscape import smooth_path
from scipy.misc import factorial
from pele.transition_states._zeroev import orthogonalize
from pele.takestep.generic import TakestepSlice

def isClose(a, b, rel_tol=1e-9, abs_tol=0.0, method='weak'):
    """
    code imported from math.isclose python 3.5
    """
    if method not in ("asymmetric", "strong", "weak", "average"):
        raise ValueError('method must be one of: "asymmetric",'
                         ' "strong", "weak", "average"')

    if rel_tol < 0.0 or abs_tol < 0.0:
        raise ValueError('error tolerances must be non-negative')

    if a == b:  # short-circuit exact equality
        return True
    # use cmath so it will work with complex or float
    if cmath.isinf(a) or cmath.isinf(b):
        # This includes the case of two infinities of opposite sign, or
        # one infinity and one finite number. Two infinities of opposite sign
        # would otherwise have an infinite relative tolerance.
        return False
    diff = abs(b - a)
    if method == "asymmetric":
        return (diff <= abs(rel_tol * b)) or (diff <= abs_tol)
    elif method == "strong":
        return (((diff <= abs(rel_tol * b)) and
                 (diff <= abs(rel_tol * a))) or
                (diff <= abs_tol))
    elif method == "weak":
        return (((diff <= abs(rel_tol * b)) or
                 (diff <= abs(rel_tol * a))) or
                (diff <= abs_tol))
    elif method == "average":
        return ((diff <= abs(rel_tol * (a + b) / 2) or
                (diff <= abs_tol)))
    else:
        raise ValueError('method must be one of:'
                         ' "asymmetric", "strong", "weak", "average"')
@jit
def compare_exact(x1, x2,
                  rel_tol=1e-9,
                  abs_tol=0.0,
                  method='weak',
                  debug=False):
    N = x1.size
    if debug:
        assert x1.size == x2.size
        assert isClose(np.dot(x1,x1), N)
        assert isClose(np.dot(x2,x2), N)
    dot = np.dot(x1, x2)
    return isClose(dot, N, rel_tol=rel_tol, abs_tol=abs_tol, method=method)

@jit
def normalize_spins(x):
    x /= (np.linalg.norm(x)/np.sqrt(len(x)))
    return x

@jit
def spin_distance_1d(x1, x2):
    x1 = normalize_spins(x1)
    x2 = normalize_spins(x2)
    return np.linalg.norm(x1 - x2)

@jit
def spin_mindist_1d(x1, x2):
    x1 = normalize_spins(x1)
    x2 = normalize_spins(x2)
    return np.linalg.norm(x1-x2), x1, x2

class UniformPSpinSPhericalRandomDisplacement(TakestepSlice):
    
    def __init__(self, nspins, stepsize=0.5):
        TakestepSlice.__init__(self, stepsize=stepsize)
        self.nspins = nspins

    def takeStep(self, coords, **kwargs):
        assert len(coords) == self.nspins
        coords[self.srange] += np.random.uniform(low=-self.stepsize, high=self.stepsize, size=coords[self.srange].shape)
        coords[self.srange] /= np.linalg.norm(coords[self.srange])/np.sqrt(self.nspins)

class MeanFieldPSpinSphericalSystem(BaseSystem):
    def __init__(self, nspins, p=3):
        BaseSystem.__init__(self)
        self.nspins = nspins
        self.p = p
        self.interactions = self.get_interactions(self.nspins, self.p)
        self.pot = self.get_potential()
        self.zerov = None
        self.setup_params(self.params)

    def setup_params(self, params):
#        params.takestep.stepsize = np.pi# / 2.
        params.takestep.verbose = True
#        nebparams.interpolator = interpolate_spins
        nebparams = params.double_ended_connect.local_connect_params.NEBparams
        nebparams.image_density = 1.2
        nebparams.iter_density = 50.
        nebparams.reinterpolate = 50
        nebparams.adaptive_nimages = True
        nebparams.adaptive_niter = True #True
        nebparams.adjustk_freq = 10
        nebparams.k = 2000
#        nebparams.distance = spin3d_distance
        params.structural_quench_params.tol = 1e-6
        params.database.overwrite_properties = False
        
        params.basinhopping.insert_rejected = True
        params.basinhopping.temperature = 10000
        
        tsparams = params.double_ended_connect.local_connect_params.tsSearchParams
        tsparams.hessian_diagonalization=True
        print "basinhopping t", params.basinhopping
        

    def get_system_properties(self):
        return dict(potential="PSpinSPherical model",
                    nspins=self.nspins,
                    p=self.p,
                    interactions=self.interactions,
                    )
        
    def get_interactions(self, nspins, p):
        assert p==3, "the interaction matrix setup at the moment requires that p==3"
        interactions = np.empty([nspins for i in xrange(p)])
        for i in xrange(nspins):
            for j in xrange(i, nspins):
                for k in xrange(j, nspins):
                    w = np.random.normal(0, np.sqrt(factorial(p)))
                    interactions[i][j][k] = w
                    interactions[k][i][j] = w
                    interactions[k][j][i] = w
                    interactions[j][k][i] = w
                    interactions[i][k][j] = w
                    interactions[j][i][k] = w
        return interactions.flatten()

    def get_potential(self, tol=1e-6):
        try:
            return self.pot
        except AttributeError:
            self.pot = MeanFieldPSpinSpherical(self.interactions, self.nspins, self.p, tol=tol)
            return self.pot

#    def _find_zero_modes(self, coords):
#        hess = self.pot.getHessian(coords)
#        wlist, vlist = np.linalg.eig(hess)
#        zerov = [v for (w,v) in  zip(wlist,vlist) if abs(w)<1e-7]
#        assert len(zerov) == 1
##        if len(zerov) > 0:
##            print "found {} - zero eigenvalue".format(len(zerov))
##            print zerov[0][:3]
#        return zerov
#
#    def _orthog_to_zero(self, v, coords):
#        if self.zerov is None:
#            self.zerov = self._find_zero_modes(coords)
#        return orthogonalize(v, self.zerov)
#    
    def get_orthogonalize_to_zero_eigenvectors(self):
        return None
        #return self._orthog_to_zero
    
    def get_metric_tensor(self, coords):
        return None
    
    def get_nzero_modes(self):
        return 1

    def get_pgorder(self, coords):
        return 1
    
    def get_mindist(self):
        return spin_mindist_1d

    def get_compare_exact(self):
        """
        are they the same minima?
        """
        return lambda x1, x2 : compare_exact(x1, x2, rel_tol=1e-7, debug=True)

    def smooth_path(self, path, **kwargs):
        mindist = self.get_mindist()
        return smooth_path(path, mindist, **kwargs)

    def get_random_configuration(self):
        coords = np.random.normal(0, 1, self.nspins)
        return normalize_spins(coords)

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
        #from pele.takestep import RandomDisplacement
        kwargs = dict(self.params["takestep"].items() + kwargs.items())
        try:
            stepsize = kwargs.pop("stepsize")
        except KeyError:
            stepsize = np.sqrt(self.nspins)/2
        #takeStep = RandomDisplacement(stepsize=stepsize)
        takeStep = UniformPSpinSPhericalRandomDisplacement(self.nspins, stepsize=stepsize)
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
        
    
    
