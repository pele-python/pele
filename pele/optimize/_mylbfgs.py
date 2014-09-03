import numpy as np
import logging
from collections import namedtuple

#from bfgs import lineSearch, BFGS
from pele.optimize import LBFGS
from mylbfgs_updatestep import mylbfgs_updatestep

__all__ = ["MYLBFGS"]

_logger = logging.getLogger("pele.optimize")

class MYLBFGS(LBFGS):
    """
    minimize a function using the LBFGS routine
    
    this class inherits everything from the pythonic LBFGS
    except the implementation of the LBFGS algorithm (determining an appropriate
    step size and direction). This is reimplemented 
    using the Fortran code from GMIN.
    
    Parameters
    ----------
    see base class LBFGS
    
    See Also
    --------
    LBFGS : base class
    """
    def __init__(self, X, pot, **lbfgs_py_kwargs):
        super(MYLBFGS, self).__init__(X, pot, **lbfgs_py_kwargs)
        
        
        N = self.N
        M = self.M
        
        #in fortran mylbfgs H0 is a vector of length N with all the elements the same
        self.H0vec = np.ones(N) * self.H0 # initial guess for the hessian

        self.W = np.zeros(N * (2 * M + 1) + 2 * M) # mylbfgs working space
        self._iter = 0
        self._point = 0
        
    
    def getStep(self, X, G):
        """
        use the Fortran MYLBFGS code from GMIN to update the LBFGS memory and get the step size and direction.
        """
        self.X = X
        self.G = G
        #save the position and gradient change
        if self._iter > 0:
            N = self.N
            M = self.M
            
            #
            # need write the change in position and the change in gradient
            # to W.  This should be done more elegantly
            #
            ISPT= N + 2*M     # index for storage of search steps
            IYPT= ISPT + N*M  # index for storage of gradient differences

            NPT = N*((self._point + M - 1) % M)  
            #s = X - self.Xold
            #y = G - self.Gold
            #print "YS YY py", np.dot( y, s ), np.dot( y,y ), ISPT+NPT
            self.W[ISPT+NPT : ISPT+NPT +N] = self.dXold
            self.W[IYPT+NPT : IYPT+NPT +N] = self.dGold    
#        self.Xold = X.copy()
#        self.Gold = G.copy()

        
        #print self._iter, self._point
        self.stp = mylbfgs_updatestep(self._iter, self.M, G, self.W, self.H0vec, self._point)
        
        #print "stp", np.linalg.norm(self.stp), self.stp
        #print "G", self.G
        #print "overlap", np.dot(self.stp, self.G)
        #print "H0", self.H0
        self.H0 = self.H0vec[0]
        self._iter += 1
        self._point = self._iter % self.M
        
        return self.stp

    def get_state(self):
        State = namedtuple("State", "W dXold dGold iter point H0 have_dXold")
        state = State(W=self.W.copy(), dXold=self.dXold.copy(), dGold=self.dGold.copy(),
                      iter=self._iter, point=self._point, H0=self.H0vec[0],
                      have_dXold=self._have_dXold
                      )
        return state
    
    def set_state(self, state):
        self.W = state.W
        self.dXold = state.dXold
        self.dGold = state.dGold
        self._have_dXold = state.have_dXold
        self._iter = state.iter
        self._point = state.point
        self.H0 = state.H0
        self.H0vec = np.ones(self.N) * self.H0
        assert self.dXold.shape == (self.N,)
        assert self.dGold.shape == (self.N,)
        assert self.W.size == (self.N * (2 * self.M + 1) + 2 * self.M)
    
    def reset(self):
        LBFGS.reset(self)
        self._iter = 0
        



#
#only testing stuff below here
# 
    

def test(pot, natoms = 100, iprint=-1): # pragma: no cover
    #import bfgs
    
    
    #X = bfgs.getInitialCoords(natoms, pot)
    #X += np.random.uniform(-1,1,[3*natoms]) * 0.3
    X = np.random.uniform(-1,1,[natoms*3])*(1.*natoms)**(1./3)*1.
    
    runtest(X, pot, natoms, iprint)

def runtest(X, pot, natoms = 100, iprint=-1): # pragma: no cover
    from _lbfgs_py import PrintEvent
    tol = 1e-5
    maxstep = 0.005

    Xinit = np.copy(X)
    e, g = pot.getEnergyGradient(X)
    print "energy", e
    
    lbfgs = LBFGS(X, pot, maxstep = 0.1, nsteps=10000, tol=tol,
                  iprint=iprint, H0=2.)
    printevent = PrintEvent( "debugout.xyz")
    lbfgs.attachEvent(printevent)
    
    ret = lbfgs.run()
    print ret
    
    print ""
    print "now do the same with scipy lbfgs"
    from pele.optimize import lbfgs_scipy as quench
    ret = quench(Xinit, pot, tol = tol)
    print ret
    #print ret[1], ret[2], ret[3]    
    
    if False:
        print "now do the same with scipy bfgs"
        from pele.optimize import bfgs as oldbfgs
        ret = oldbfgs(Xinit, pot, tol = tol)
        print ret
    
    if False:
        print "now do the same with gradient + linesearch"
        import _bfgs
        gpl = _bfgs.GradientPlusLinesearch(Xinit, pot, maxstep = 0.1)  
        ret = gpl.run(1000, tol = 1e-6)
        print ret
            
    if False:
        print "calling from wrapper function"
        from pele.optimize import lbfgs_py
        ret = lbfgs_py(Xinit, pot, tol = tol)
        print ret
        
    if True:
        print ""
        print "now do the same with lbfgs_py"
        from pele.optimize import lbfgs_py
        ret = lbfgs_py(Xinit, pot, tol = tol)
        print ret



    try:
        import pele.utils.pymolwrapper as pym
        pym.start()
        for n, coords in enumerate(printevent.coordslist):
            coords=coords.reshape(natoms, 3)
            pym.draw_spheres(coords, "A", n)
    except ImportError:
        print "error loading pymol"

def test2(): # pragma: no cover
    #from pele.potentials.lj import LJ as Pot
    from pele.potentials.ATLJ import ATLJ as Pot
    pot = Pot()

    test(pot, natoms=3, iprint=-1)
    exit(1)
    
    coords = np.loadtxt("coords")
    print coords.size
    coords = np.reshape(coords, coords.size)
    print coords
    runtest(coords, pot, natoms=3, iprint=1)
    
if __name__ == "__main__":
    test2()
