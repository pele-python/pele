import numpy as np
from bfgs import lineSearch, BFGS
import lbfgs_py
from mylbfgs_updatestep import mylbfgs_updatestep

class LBFGS(lbfgs_py.LBFGS):
    def __init__(self, X, pot, **lbfgs_py_kwargs):
        """
        this class inherits everything from the pythonic LBFGS
        except the implementation of the LBFGS (determining an appropriate
        step size and direction). This is reimplemented 
        using the Fortran code from GMIN.
        """
        lbfgs_py.LBFGS.__init__(self, X, pot, **lbfgs_py_kwargs)
        
        
        N = self.N
        M = self.M
        
        #js850> why is H0 length N here and length M in lbfgs_py?
        self.H0 = np.ones(N) * 1. #initial guess for the hessian

        self.W = np.zeros(N*(2*M+1)+2*M) #mylbfgs working space
        self.iter = 0
        self.point = 0
        
    
    def getStep(self, X, G):
        """
        use the Fortran LBFGS code from GMIN to get the step
        size and direction.
        """
        self.X = X
        self.G = G
        #save the position and gradient change
        if self.iter > 0:
            N = self.N
            M = self.M
            
            #
            # need write the change in position and the change in gradient
            # to W.  This should be done more elegantly
            #
            ISPT= N + 2*M     # index for storage of search steps
            IYPT= ISPT + N*M  # index for storage of gradient differences

            NPT = N*((self.point + M - 1) % M)  
            s = X - self.Xold
            y = G - self.Gold
            #print "YS YY py", np.dot( y, s ), np.dot( y,y ), ISPT+NPT
            self.W[ISPT+NPT : ISPT+NPT +N] = X - self.Xold
            self.W[IYPT+NPT : IYPT+NPT +N] = G - self.Gold    
        self.Xold = X.copy()
        self.Gold = G.copy()

        
        #print self.iter, self.point
        self.stp = mylbfgs_updatestep(self.iter, self.M, G, self.W, self.H0, self.point, [self.N])
        
        #print "stp", np.linalg.norm(self.stp), self.stp
        #print "G", self.G
        #print "overlap", np.dot(self.stp, self.G)
        #print "H0", self.H0
        self.iter += 1
        self.point = self.iter % self.M
        
        return self.stp

        
    

def test(pot, natoms = 100, iprint=-1):
    import bfgs
    
    
    #X = bfgs.getInitialCoords(natoms, pot)
    #X += np.random.uniform(-1,1,[3*natoms]) * 0.3
    X = np.random.uniform(-1,1,[natoms*3])*(1.*natoms)**(1./3)*1.
    
    runtest(X, pot, natoms, iprint)

def runtest(X, pot, natoms = 100, iprint=-1):
    from lbfgs_py import PrintEvent
    tol = 1e-5
    maxstep = 0.005

    Xinit = np.copy(X)
    e, g = pot.getEnergyGradient(X)
    print "energy", e
    
    lbfgs = LBFGS(X, pot, maxstep = 0.1, nsteps=10000, tol=tol,
                  iprint=iprint)
    printevent = PrintEvent( "debugout.xyz")
    lbfgs.attachEvent(printevent)
    
    ret = lbfgs.run()
    print "done", ret[1], ret[2], ret[3], ret[5]
    
    print ""
    print "now do the same with scipy lbfgs"
    from pygmin.optimize.quench import quench
    ret = quench(Xinit, pot.getEnergyGradient, tol = tol)
    print ret[1], ret[2], ret[3]    
    
    if False:
        print "now do the same with scipy bfgs"
        from pygmin.optimize.quench import bfgs as oldbfgs
        ret = oldbfgs(Xinit, pot.getEnergyGradient, tol = tol)
        print ret[1], ret[2], ret[3]    
    
    if False:
        print "now do the same with gradient + linesearch"
        import bfgs
        gpl = bfgs.GradientPlusLinesearch(Xinit, pot, maxstep = 0.1)  
        ret = gpl.run(1000, tol = 1e-6)
        print ret[1], ret[2], ret[3]    
            
    if False:
        print "calling from wrapper function"
        from pygmin.optimize.quench import lbfgs_py
        ret = lbfgs_py(Xinit, pot.getEnergyGradient, tol = tol)
        print ret[1], ret[2], ret[3]    
        
    if True:
        print ""
        print "now do the same with lbfgs_py"
        from pygmin.optimize.quench import lbfgs_py
        ret = lbfgs_py(Xinit, pot.getEnergyGradient, tol = tol)
        print ret[1], ret[2], ret[3]    



    if True:
        import pygmin.utils.pymolwrapper as pym
        pym.start()
        for n, coords in enumerate(printevent.coordslist):
            coords=coords.reshape(natoms, 3)
            pym.draw_spheres(coords, "A", n)

        
if __name__ == "__main__":
    #from pygmin.potentials.lj import LJ as Pot
    from pygmin.potentials.ATLJ import ATLJ as Pot
    pot = Pot()

    test(pot, natoms=3, iprint=-1)
    exit(1)
    
    coords = np.loadtxt("coords")
    print coords.size
    coords = np.reshape(coords, coords.size)
    print coords
    runtest(coords, pot, natoms=3, iprint=1)
    












