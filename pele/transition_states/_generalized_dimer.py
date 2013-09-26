import numpy as np

from pele.transition_states import findLowestEigenVector
from pele.optimize import mylbfgs, Result

def find_TS_generalized_dimer(coords, potential, eigenvec0=None, minimizer=mylbfgs,
                              dimer_kwargs=None, minimizer_kwargs=None):
    if dimer_kwargs is None: dimer_kwargs = {}
    if minimizer_kwargs is None: minimizer_kwargs = {}
    if eigenvec0 is None:
        eigenvec0 = np.random.rand(coords.size)
    eigenvec0 /= np.linalg.norm(eigenvec0)
    assert coords.shape == eigenvec0.shape

    dimer = GeneralizedDimer(potential, eigenvec0, **dimer_kwargs)
    
    qres = minimizer(coords, dimer, **minimizer_kwargs)
    
    res = Result()
    res.eigenval = dimer.eigenval
    res.eigenvec = dimer.eigenvec 
    res.coords = qres.coords
    res.energy = dimer.energy
    res.grad = dimer.true_gradient
    res.rms = qres.rms
    res.nfev = dimer.nfev
    res.nsteps = qres.nsteps
    
    return res
    

class GeneralizedDimer(object):
    """
    """
    def __init__(self, potential, eigenvec0, find_lowest_eigenvector=findLowestEigenVector, leig_tol=1e-4,
                 n_translational_steps=5, n_rotational_steps=20, H0_evec=None):
        self.potential = potential
        self.find_lowest_eigenvector = find_lowest_eigenvector
        self.eigenvec = eigenvec0 / np.linalg.norm(eigenvec0)
        self.leig_tol = leig_tol
    
        self.n_translational_steps = n_translational_steps
        self.n_rotational_steps = n_rotational_steps
        self.iter_number = 0
        
        self._H0 = H0_evec
        self.nfev = 0
    
    def getEnergyGradientInverted(self, x):
        e, g = self.potential.getEnergyGradient(x)
        self.energy = e
        self.true_gradient = g.copy()
        g -= 2. * np.dot(g, self.eigenvec) * self.eigenvec
        self.nfev += 1
        return e, g

    def update_eigenvec(self, x):
        ret = self.find_lowest_eigenvector(x, self.potential, eigenvec0=self.eigenvec, H0=self._H0, 
                                           dx=1e-3, 
                                           nsteps=self.n_rotational_steps, tol=self.leig_tol)
        self._H0 = ret.H0
        self.eigenvec = ret.eigenvec
        self.eigenval = ret.eigenval
        self.nfev += ret.nfev
    
    def getEnergyGradient(self, x):
        """this is the main loop of the program.  it will be called by the optimizer
        
        """
        if self.iter_number % self.n_translational_steps == 0:
            self.update_eigenvec(x)
        
        self.iter_number += 1
        e, g = self.getEnergyGradientInverted(x)
        return 0., g
        
#
# testing only below here
#
        
class PotWrapper(object):
    def __init__(self, pot):
        self.pot = pot
        self.nfev = 0
    
    def getEnergyGradient(self, x):
        self.nfev += 1
        return self.pot.getEnergyGradient(x)
        
def compare_HEF(x0, evec0, system, **kwargs):
    from pele.transition_states import findTransitionState
    pot = PotWrapper(system.get_potential())
    ret = findTransitionState(x0, pot, eigenvec0=evec0, orthogZeroEigs=None, **kwargs)
    print ret.eigenval
    print pot.nfev
    print ret.rms


def get_x0():
    from pele.systems import LJCluster
    natoms = 31
    system = LJCluster(natoms)
    db = system.create_database()
    bh = system.get_basinhopping(db, outstream=None)
    while db.number_of_minima() < 2:
        bh.run(1)
    
    mindist = system.get_mindist()
    m1, m2 = db.minima()[:4]
    d, x1, x2 = mindist(m1.coords, m2.coords)
    
    x0 = (x1 + x2) / 2
    evec0 = x2 - x1
    
    return system, x0, evec0


def test():
    system, x0, evec0 = get_x0() 
    
    ret = find_TS_generalized_dimer(x0.copy(), system.get_potential(), 
                                    eigenvec0=evec0, 
                                    dimer_kwargs=dict(n_translational_steps=5,
                                                      n_rotational_steps=20), 
                                    minimizer_kwargs=dict(iprint=1, tol=4e-5, nsteps=2000)
                                    )
    
    print "eigenvalue", ret.eigenval
    print "function evaluations", ret.nfev
    print "rms", ret.rms
    
    print "\n\nnow with hybrid eigenvector following"
    compare_HEF(x0, evec0, system, nsteps_tangent1=5, nsteps_tangent2=5, 
                lowestEigenvectorQuenchParams={"nsteps":20, "tol":1e-4}, iprint=10)
    
    
        

if __name__ == "__main__":
    test()
