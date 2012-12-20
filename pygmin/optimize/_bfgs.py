import numpy as np
#from itertools import xrange
from pygmin.potentials.potential import potential as basepotential

class LineSearchPot(basepotential):
    """
    project the N dimensional potential into a one dimensional 
    problem along the line defined by X + a*V
    """
    def __init__(self, X, V, pot):
        self.X = X.copy()        
        self.V = V.copy()
        self.pot = pot
        self.Vnorm = np.linalg.norm(self.V)
        self.V /= self.Vnorm
    def getEnergy(self, a):
        return self.pot.getEnergy( self.X + a * self.V)
    def getEnegyGradient(self, a):
        e, g = self.pot.getEnergyGradient( self.X + a * self.V )
        #project g onto V
        g1 = np.dot(self.V, g)
        return e, g1

def lineSearch(X, V, pot, aguess = 0.1, tol = 1e-3):
    """
    minimize the potential along the line defined by X + a * V
    
    i'm sure there's a better way to do this
    """
    tol /= np.sqrt(len(X)) #normalize the tolerance because we're reducing dimensionality
    ls = LineSearchPot(X, V, pot)
    
    einit = pot.getEnergy(X)
    #from optimize.quench import steepest_descent as quench
    from pygmin.optimize import lbfgs_scipy as quench
    a = np.zeros(1) * aguess
    ret = quench(a, ls.getEnergyGradient, tol=tol)
    a = ret[0][0]
    e = ret[1]
    funcalls = ret[3] 
    
    #print "    linesearch step size", a, e, einit, e - einit, funcalls
    return a, e, funcalls

def lineSearchScipy(X, V, pot, aguess = 0.1, tol = 1e-3, G = None):
    """
    js850> this doesn't seem to work
    """
    import scipy.optimize
    V = V.copy() 
    V /= np.linalg.norm(V)
    #ret = scipy.optimize.line_search(pot.getEnergy, pot.getGradient, X, V, gfk = G, old_fval = E)
    ret = scipy.optimize.line_search(pot.getEnergy, pot.getGradient, X, V, gfk = G)
    a = ret[0]
    funcalls = ret[1] + ret[2]
    return a, funcalls
    

        

class BFGS:
    def __init__(self, X, pot, maxstep):
        self.maxstep = maxstep
        self.pot = pot
        self.X = X
        self.funcalls = 0
    
        self.N = len(X)
        N = self.N
        
        self.B = np.zeros([N,N]) #the approximate hessian
        #self.Bprev = np.copy(self.B)
        self.Gprev = np.zeros(N)
        self.k = 0
        
        
        #self.B[:] = 1. / (G[:, np.newaxis] * G[np.newaxis, :] )
        self.B[:] = np.identity(N)
                
        self.stp = np.zeros(N)
    
    def step(self, X, G):
        self.G = G
        B = self.B
        
        if self.k > 0:
            #update B
            s = self.stp
            y = G - self.Gprev
            
            B += (np.outer(y,y) / np.dot(y,s)) - \
                (np.dot( B, np.dot( np.outer(s,s), B) ) / np.dot( s, np.dot(B, s)))

        
        #solve for p in B*p = G
        self.stp = -np.linalg.solve(B, G)
        
        self.takeStep(X, self.stp)
        
            
        
        
        self.k += 1
        
        self.Gprev[:] = G[:]
        #self.Bprev[:] = B[:] 
        return X

    def takeStep(self, X, stp):
        self.stp = stp
        #we have the step direction.  now take the step
        stepsize = self.getStepSize(X, self.stp)
        snorm = np.linalg.norm(self.stp)
        self.stp *= stepsize / snorm
        smax = np.max(np.abs(self.stp))
        if True and smax > self.maxstep:
            #print "reducing step from", smax, "to", self.maxstep
            self.stp *= self.maxstep / smax
        #print "    step size", np.linalg.norm(self.stp)
        X += self.stp

    
    def getStepSize(self, X, stp):
        aguess = np.linalg.norm(stp)
        a, e, funcalls = lineSearch(X, stp, self.pot, aguess = aguess, tol = self.tol*.05)
        #a2, funcalls2 = lineSearchScipy(X, stp, self.pot, aguess = aguess, tol = self.tol*.05, G = self.G)
        #print "step size", a, a2, funcalls, funcalls2
        self.funcalls += funcalls
        return a
    
        
    
    def run(self, nsteps = 1000, tol = 1e-6, iprint = -1):
        self.tol = tol
        X = self.X
        sqrtN = np.sqrt(self.N)
        
        i = 1
        self.funcalls += 1
        e, G = self.pot.getEnergyGradient(X)
        while i < nsteps:
            X = self.step(X, G)
            e, G = self.pot.getEnergyGradient(X)
            
            rms = np.linalg.norm(G) / sqrtN

            i += 1
            self.funcalls += 1
            
            if iprint > 0:
                if i % iprint == 0:
                    print "quench step", i, e, rms, self.funcalls
      
            if rms < self.tol:
                break
        return X, e, rms, self.funcalls, G , i
    
class GradientPlusLinesearch(BFGS):
    """
    this is only for testing.  Anything using linesearch should run better than this.
    """
    def __init__(self, X, pot, maxstep):
        self.maxstep = maxstep
        self.pot = pot
        self.X = X
        self.funcalls = 0
        self.N = len(X)
    
    def step(self, X, G):
        self.G = G
        #solve for p in B*p = G
        self.stp = -G
        
        
        #we have the step direction.  now take the step
        stepsize = self.getStepSize(X, self.stp)
        snorm = np.linalg.norm(self.stp)
        self.stp *= stepsize / snorm
        smax = np.max(np.abs(self.stp))
        if True and smax > self.maxstep:
            #print "reducing step from", smax, "to", self.maxstep
            self.stp *= self.maxstep / smax
        #print "    step size", np.linalg.norm(self.stp)
        X += self.stp
        
        return X


    


def getInitialCoords(natoms, pot):
    from pygmin.basinhopping import BasinHopping
    from pygmin.takestep.displace import RandomDisplacement
    takestep = RandomDisplacement(0.3)
    X = np.random.uniform(-1,1,natoms*3)*(1.*natoms)**(1./3)*.1
    bh = BasinHopping(X, pot, takestep)
    bh.run(30)
    
    X = bh.coords
    return X



def test():
    natoms = 100
    tol = 1e-6
    
    from pygmin.potentials.lj import LJ
    pot = LJ()
    
    X = getInitialCoords(natoms, pot)
    X += np.random.uniform(-1,1,[3*natoms]) * 0.3
    
    #do some steepest descent steps so we don't start with a crazy structure
    #from optimize.quench import _steepest_descent as steepestDescent
    #ret = steepestDescent(X, pot.getEnergyGradient, iprint = 1, dx = 1e-4, nsteps = 100, gtol = 1e-3, maxstep = .5)
    #X = ret[0]
    #print X

    
    Xinit = np.copy(X)
    e, g = pot.getEnergyGradient(X)
    print "energy", e
    
    lbfgs = BFGS(X, pot, maxstep = 0.1)
    
    ret = lbfgs.run(100, tol = tol, iprint=1)
    print "done", ret[1], ret[2], ret[3]
    
    print "now do the same with scipy lbfgs"
    from pygmin.optimize import lbfgs_scipy as quench
    ret = quench(Xinit, pot.getEnergyGradient, tol = tol)
    print ret[1], ret[2], ret[3]    
    
    print "now do the same with scipy bfgs"
    from pygmin.optimize import bfgs as oldbfgs
    ret = oldbfgs(Xinit, pot.getEnergyGradient, tol = tol)
    print ret[1], ret[2], ret[3]    
    
    print "now do the same with old gradient + linesearch"
    gpl = GradientPlusLinesearch(Xinit, pot, maxstep = 0.1)  
    ret = gpl.run(100, tol = 1e-6)
    print ret[1], ret[2], ret[3]    
    


        
if __name__ == "__main__":
    test()
    












