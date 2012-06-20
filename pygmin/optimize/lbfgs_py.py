import numpy as np
from bfgs import lineSearch, BFGS

class LBFGS:
    def __init__(self, X, pot, maxstep = 0.1, maxErise = 1e-4, M=10):
        self.X = X
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
        
        self.s = np.zeros([M,N])  #position updates
        self.y = np.zeros([M,N])  #gradient updates
        self.a = np.zeros(M)  #approximation for the inverse hessian
        #self.beta = np.zeros(M) #working space
        
        self.q = np.zeros(N)  #working space
        
        self.H0 = np.ones(M) * 1. #initial guess for the hessian
        self.rho = np.zeros(M)
        self.k = 0
        
        self.s[0,:] = self.X
        self.y[0,:] = self.G
        self.rho[0] = 0. #1. / np.dot(X,G)
        
        self.stp = np.zeros(N)

        self.Xold = self.X.copy()
        self.Gold = self.G.copy()
        
        self.nfailed = 0
    
    def step(self, X, G):
        """
        see http://en.wikipedia.org/wiki/Limited-memory_BFGS
        """
        self.G = G #saved for the line search
        
        s = self.s
        y = self.y
        a = self.a
        q = self.q
        rho = self.rho
        M = self.M
        
        k = self.k
        ki = k % M #the index corresponding to k
        
        
        #we have a new X and G, save in s and y
        if k > 0:
            km1 = (k + M - 1) % M  #=k-1  cyclical
            s[km1,:] = X - self.Xold
            y[km1,:] = G - self.Gold
            
            YS = np.dot(s[km1,:], y[km1,:])
            if YS == 0.:
                print "warning: resetting YS to 1 in lbfgs", YS
                YS = 1.            
            rho[km1] = 1. / YS
            
            #update the approximation for the diagonal inverse hessian
            YY = np.dot( y[km1,:], y[km1,:] )
            if YY == 0.:
                print "warning: resetting YY to 1 in lbfgs", YY
                YY = 1.
            self.H0[ki] = YS / YY

        self.Xold[:] = X[:]
        self.Gold[:] = G[:]

        
        q[:] = G[:]
        myrange = [ i % M for i in range(max([0,k-M]), k, 1) ]
        #print "myrange", myrange, ki, k
        for i in reversed(myrange):
            a[i] = rho[i] * np.dot( s[i,:], q )
            q -= a[i] * y[i,:]
        
        #z[:] = self.H0[ki] * q[:]
        z = q #q is not used anymore, so we can use it as workspace
        z *= self.H0[ki]
        for i in myrange:
            beta = rho[i] * np.dot( y[i,:], z )
            z += s[i,:] * (a[i] - beta)
        
        self.stp[:] = -z[:]
        
        #we now have the step direction.  now take the step
        #self.takeStep(X, self.stp)
        #print "step size", np.linalg.norm(self.stp)
        
        self.k += 1
        return self.stp

    def takeStepNoLineSearch(self, X, E, G, stp):
        """
        We now have a proposed step.  This function will make sure it is 
        a good step and then take it.
        
        1) if the step is not anti-aligned with the gradient (i.e. downhill), then reverse the step
        
        2) if the step is larger than maxstep, then rescale the step
        
        3) calculate the energy and gradient of the new position
        
        4) if the step increases the energy by more than maxErise, 
            then reduce the step size and go to 3)
                
        *The below is not implemented yet.  It's on the TODO list
        
        6) if failures occur too many times, restart the quench process from the current configuration
        
        7) if we're still failing then abort
        """
        f = 1.
        X0 = X.copy()
        G0 = G.copy()
        E0 = E
        maxErise = self.maxErise
        
        if np.dot(G, stp) > 0:
            #print "overlap was negative, reversing step direction"
            stp = -stp
            self.stp = stp
        
        stepsize = np.linalg.norm(stp)
        
        if f*stepsize > self.maxstep:
            f = self.maxstep / stepsize
        
        nincrease = 0
        while True:
            X = X0 + f * stp
            E, G = self.pot.getEnergyGradient(X)
            self.funcalls += 1
            
            if (E - E0) <= maxErise:
                break
            else:
                #print "warning: energy increased, trying a smaller step", E, E0, f*stepsize
                f /= 10.
                nincrease += 1

        if nincrease <= 1:
            self.nfailed = 0
        else:
            self.nfailed += 1
            if False and self.nfailed > 3:
                print "resetting H0"
                print self.H0
                self.reset()
        
        if False and self.k <= 1:
            print G0
            print stp
            print G
            
        
        return X, E, G
    
    def reset(self):
        self.H0[:] = 1.
        self.k = 0
    
    def attachEvent(self, event):
        self.events.append(event)
                
    def run(self, nsteps = 1000, tol = 1e-6, iprint = -1):
        self.tol = tol
        X = self.X
        sqrtN = np.sqrt(self.N)
        
        
        i = 1
        self.funcalls += 1
        e, G = self.pot.getEnergyGradient(X)
        while i < nsteps:
            stp = self.step(X, G)
            
            X, e, G = self.takeStepNoLineSearch(X, e, G, stp)
            #e, G = self.pot.getEnergyGradient(X)
            
            rms = np.linalg.norm(G) / sqrtN

            i += 1
            
            if iprint > 0:
                if i % iprint == 0:
                    print "quench step", i, e, rms, self.funcalls
            for event in self.events:
                event( coords=X, energy=e, rms=rms )
      
            if rms < self.tol:
                break
        return X, e, rms, self.funcalls, G , i
   

class PrintEvent:
    def __init__(self, fname):
        self.fout = open(fname, "w")
        self.coordslist = []

    def __call__(self, coords, **kwargs):
        from pygmin.printing.print_atoms_xyz import printAtomsXYZ as printxyz 
        printxyz(self.fout, coords)
        self.coordslist.append( coords.copy() )
        
    

def test(pot, natoms = 100, iprint=-1):
    import bfgs
    
    
    #X = bfgs.getInitialCoords(natoms, pot)
    #X += np.random.uniform(-1,1,[3*natoms]) * 0.3
    X = np.random.uniform(-1,1,[natoms*3])*(1.*natoms)**(1./3)*.5
    
    runtest(X, pot, natoms, iprint)

def runtest(X, pot, natoms = 100, iprint=-1):
    tol = 1e-5

    Xinit = np.copy(X)
    e, g = pot.getEnergyGradient(X)
    print "energy", e
    
    lbfgs = LBFGS(X, pot, maxstep = 0.1)
    printevent = PrintEvent( "debugout.xyz")
    lbfgs.attachEvent(printevent)
    
    ret = lbfgs.run(10000, tol = tol, iprint=iprint)
    print "done", ret[1], ret[2], ret[3], ret[5]
    
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
        import pygmin.utils.pymolwrapper as pym
        pym.start()
        for n, coords in enumerate(printevent.coordslist):
            coords=coords.reshape(natoms, 3)
            pym.draw_spheres(coords, "A", n)

        
if __name__ == "__main__":
    from pygmin.potentials.lj import LJ
    from pygmin.potentials.ATLJ import ATLJ
    pot = ATLJ()

    #test(pot, natoms=3, iprint=1)
    
    coords = np.loadtxt("coords")
    print coords.size
    coords = np.reshape(coords, coords.size)
    print coords
    runtest(coords, pot, natoms=3, iprint=1)
    












