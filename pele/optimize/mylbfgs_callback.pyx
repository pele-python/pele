import numpy as np
from bfgs import lineSearch, BFGS
import lbfgs_py
import mylbfgs_fort

class LBFGS(object):
    def __init__(self, X, pot, maxstep = 0.1, maxErise = 1e-4, M=10):
        self.X = X.copy()
        self.pot = pot
        e, self.G = self.pot.getEnergyGradient(self.X)
        self.funcalls = 1
        self.maxstep = maxstep
        self.maxErise = maxErise
        self.events = []
    
        self.N = len(X)
        self.M = M 
        N = self.N
        M = self.M
        
        #self.beta = np.zeros(M) #working space
        
    def getEnergyGradientCallback(self, coords):
        e, grad = self.pot.getEnergyGradient(coords)
        #print "in getEnergyGradientCallback", e
        return grad, e
    
    
    def run(self, nsteps = 1000, tol = 1e-6, iprint = -1):
        ret = mylbfgs_fort.mylbfgs(self.M, self.X, tol, nsteps, self.maxstep, \
                self.maxErise, self.getEnergyGradientCallback)
        
        X, mflag, e, itdone, funcalls, G = ret
        
        rms = np.linalg.norm(G) / np.sqrt(self.N)
        return X, e, rms, funcalls, G , itdone




        
    

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
    
    lbfgs = LBFGS(X, pot, maxstep = 0.1)
    printevent = PrintEvent( "debugout.xyz")
    #lbfgs.attachEvent(printevent)
    
    ret = lbfgs.run(10000, tol = tol, iprint=iprint)
    print "done", ret[1], ret[2], ret[3], ret[5]
    
    print ""
    print "now do the same with scipy lbfgs"
    from pygmin.optimize import lbfgs_scipy as quench
    ret = quench(Xinit, pot.getEnergyGradient, tol = tol)
    print ret[1], ret[2], ret[3]    
    
    if False:
        print "now do the same with scipy bfgs"
        from pygmin.optimize import bfgs as oldbfgs
        ret = oldbfgs(Xinit, pot.getEnergyGradient, tol = tol)
        print ret[1], ret[2], ret[3]    
    
    if False:
        print "now do the same with gradient + linesearch"
        import bfgs
        gpl = bfgs.GradientPlusLinesearch(Xinit, pot, maxstep = 0.1)  
        ret = gpl.run(1000, tol = 1e-6)
        print ret[1], ret[2], ret[3]    
            
    if True:
        print ""
        print "calling from wrapper function"
        from pygmin.optimize import mylbfgs as quench
        ret = quench(Xinit, pot.getEnergyGradient, tol = tol)
        print ret[1], ret[2], ret[3]    
        
    if True:
        print ""
        print "now do the same with lbfgs_py"
        from pygmin.optimize import lbfgs_py
        ret = lbfgs_py(Xinit, pot.getEnergyGradient, tol = tol)
        print ret[1], ret[2], ret[3]    



    if False:
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
    












