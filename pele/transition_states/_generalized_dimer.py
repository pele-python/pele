import numpy as np

from pele.transition_states import findLowestEigenVector

class GeneralizedDimer(object):
    """ 
    """
    def __init__(self, potential, find_lowest_eigenvector, eigenvec0, leig_tol=1e-4,
                 n_translational_steps=5, n_rotational_steps=20):
        self.potential = potential
        self.find_lowest_eigenvector = find_lowest_eigenvector
        self.eigenvec = eigenvec0 / np.linalg.norm(eigenvec0)
        self.leig_tol = leig_tol
    
        self.n_translational_steps = n_translational_steps
        self.n_rotational_steps = n_rotational_steps
        self.n_iter = 0
        
        self._H0 = None
        self.nfev = 0
    
    def getEnergyGradientInverted(self, x):
        e, g = self.potential.getEnergyGradient(x)
        g -= 2. * np.dot(g, self.eigenvec) * self.eigenvec
        self.energy = e
        self.gradient = g
        self.nfev += 1
        return e, g

    
    def getEnergyGradient(self, x):
        """this is the main loop of the program.  it will be called by the optimizer
        
        """
        if self.n_iter % self.n_translational_steps == 0:
            ret = self.find_lowest_eigenvector(x, self.potential, eigenvec0=self.eigenvec, H0=self._H0, 
                                               orthogZeroEigs=0, dx=1e-3, 
                                               nsteps=self.n_rotational_steps, tol=self.leig_tol)
            self._H0 = ret.H0
            self.eigenvec = ret.eigenvec
            self.eigenval = ret.eigenval
            self.nfev += ret.nfev
            
            print self.n_iter, "eval", self.eigenval

        
        self.n_iter += 1
        e, g = self.getEnergyGradientInverted(x)
        return 0., g
        
        
        
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
    from pele.optimize import mylbfgs
    system, x0, evec0 = get_x0() 
    
    tssearch = GeneralizedDimer(system.get_potential(), 
                                findLowestEigenVector, evec0)
    
    
    ret = mylbfgs(x0.copy(), tssearch, iprint=1, tol=4e-5, nsteps=2000)
    
    print "eigenvalue", tssearch.eigenval
    print "function evaluations", tssearch.nfev
    print "rms", ret.rms
    
    print "\n\nnow with hybrid eigenvector following"
    compare_HEF(x0, evec0, system, nsteps_tangent1=5, nsteps_tangent2=5, 
                lowestEigenvectorQuenchParams={"nsteps":20, "tol":1e-4})
    
    
        

if __name__ == "__main__":
    test()
