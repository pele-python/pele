import numpy as np

from pele.transition_states import FindLowestEigenVector, analyticalLowestEigenvalue
from pele.optimize import MYLBFGS, Result
from pele.utils import rotations

def _run_dimer(dimer, translator, rotator):
    # find the lowest eigenvector for the first time
    ret = rotator.run(100)
    dimer.update_eigenvec(ret.eigenvec, ret.eigenval)

    print "trans rot", dimer.n_translational_steps, dimer.n_rotational_steps
    while True:
        # translate the dimer
        print "translating dimer"
        for i in xrange(dimer.n_translational_steps):
            if translator.stop_criterion_satisfied() or translator.iter_number >= translator.nsteps:
                return
            translator.one_iteration()
        
        # update the eigenvector (rotate the dimer)
        print "rotating dimer"
        mret = translator.get_result()
        rotator.update_coords(mret.coords)#, energy=dimer.energy, gradient=dimer.true_gradient)
        ret = rotator.run(dimer.n_rotational_steps)
        dimer.update_eigenvec(ret.eigenvec, ret.eigenval)


def find_TS_generalized_dimer(coords, potential, eigenvec0=None, minimizer_class=MYLBFGS,
                              leig_kwargs=None,
                              dimer_kwargs=None, minimizer_kwargs=None):
    # check the keyword dictionaries
    if dimer_kwargs is None: dimer_kwargs = {}
    if minimizer_kwargs is None: minimizer_kwargs = {}
    if leig_kwargs is None: leig_kwargs = {}
    if eigenvec0 is None:
        eigenvec0 = rotations.vec_random_ndim(coords.shape)
    eigenvec0 /= np.linalg.norm(eigenvec0)
    assert coords.shape == eigenvec0.shape
    
    dimer_kwargs["auto_rotate"] = False

    # set up the object that will maintain the rotation of the dimer
    rotator = FindLowestEigenVector(coords, potential, eigenvec0=eigenvec0, orthogZeroEigs=0, **leig_kwargs)

    # set up the dimer potential
    dimer = GeneralizedDimer(potential, eigenvec0, **dimer_kwargs)
    
    # set up the optimizer that will translate the dimer
    translator = minimizer_class(coords, dimer, **minimizer_kwargs)
    
    # optimize the dimer
    _run_dimer(dimer, translator, rotator)
    
    qres = translator.get_result()
    
    res = Result()
    res.eigenval = dimer.eigenval
    res.eigenvec = dimer.eigenvec 
    res.coords = qres.coords
    res.energy = dimer.energy
    res.grad = dimer.true_gradient
    res.rms = qres.rms
    res.nfev = dimer.nfev
    res.nsteps = qres.nsteps
    res.success = qres.success
    
    if res.eigenval > 0:
        res.success = False
    
    return res
    

class GeneralizedDimer(object):
    """
    """
    def __init__(self, potential, eigenvec0,
                 n_translational_steps=5, n_rotational_steps=20, leig_H0=None,
                 leig_kwargs=None, auto_rotate=True,
                 ):
        self.potential = potential
        self.eigenvec = eigenvec0 / np.linalg.norm(eigenvec0)
    
        self.n_translational_steps = n_translational_steps
        self.n_rotational_steps = n_rotational_steps
        self.iter_number = 0
        self.auto_rotate = auto_rotate
        
        self._H0 = leig_H0
        self.nfev = 0
        
        self._leig_minimizer_state = None
        self.leig_kwargs = leig_kwargs
        if self.leig_kwargs is None:
            self.leig_kwargs = dict()
    
    def getEnergyGradientInverted(self, x):
        e, g = self.potential.getEnergyGradient(x)
        self.energy = e
        self.true_gradient = g.copy()
        g -= 2. * np.dot(g, self.eigenvec) * self.eigenvec
        self.nfev += 1
        return e, g

    def update_eigenvec_analytical(self, x, **kwargs):
        self.eigenval, self.eigenvec = analyticalLowestEigenvalue(x, self.potential)
        self.nfev += 1

#    def update_eigenvec(self, x, n_rotational_steps=None):
#        if n_rotational_steps is None:
#            n_rotational_steps = self.n_rotational_steps
#        ret = self.find_lowest_eigenvector(x, self.potential, eigenvec0=self.eigenvec, H0=self._H0, 
#                                           minimizer_state=self._leig_minimizer_state,
#                                           nsteps=n_rotational_steps, 
#                                           **self.leig_kwargs)
#        self._H0 = ret.H0
#        self._leig_minimizer_state = ret.minimizer_state
#        self.eigenvec = ret.eigenvec
#        self.eigenval = ret.eigenval
#        self.nfev += ret.nfev
    
    def update_eigenvec(self, eigenvec, eigenval):
        self.eigenvec = eigenvec.copy()
        self.eigenval = eigenval
    
    def getEnergyGradient(self, x):
        """this is the main loop of the program.  it will be called by the optimizer
        
        """
#        if self.iter_number % self.n_translational_steps == 0 and self.auto_rotate:
#            self.update_eigenvec(x)
        
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
