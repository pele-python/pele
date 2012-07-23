import numpy as np
from pygmin.potentials.potential import potential as basepot
from orthogopt import orthogopt

class LowestEigPot(basepot):
    """
    this is a potential wrapper designed to use optimization to find the eigenvector
    which corresponds to the lowest eigenvalue
    
    here the energy corresponds to the eigenvalue, and the coordinates to be optimized is the eigenvector
    """
    def __init__(self, coords, pot):
        """
        coords is the point in space where we want to compute the lowest eigenvector
        
        pot is the potential of the system.  
            i.e. pot.getEnergyGradient(coords) gives the energy and gradient
        """
        self.coords = np.copy(coords)
        self.pot = pot
        self.E, self.G = self.pot.getEnergyGradient(self.coords)
        
        self.diff = 1e-3
    
    
    def getEnergyGradient(self, vec_in):
        """
        vec is a guess for the lowest eigenvector.  It should be normalized
        """
        vecl = 1.
        vec = vec_in / np.linalg.norm(vec_in)
        vec = orthogopt(vec, self.coords, True)
        coordsnew = self.coords - self.diff * vec
        Eminus, Gminus = self.pot.getEnergyGradient(coordsnew)
        
        coordsnew = self.coords + self.diff * vec
        Eplus, Gplus = self.pot.getEnergyGradient(coordsnew)
        
        #diag = (Eplus + Eminus -2.0 * self.E) / (self.diff**2, vecl)
        
        diag2 = np.sum((Gplus - Gminus) * vec) / (2.0 * self.diff)
        
        """
        DIAG3=2*(DIAG-DIAG2/2)
        C  Although DIAG3 is a more accurate estimate of the diagonal second derivative, it
        C  cannot be differentiated analytically.
        """
        
        #GL(J1)=(GRAD1(J1)-GRAD2(J1))/(ZETA*VECL**2)-2.0D0*DIAG2*LOCALV(J1)/VECL**2
        grad = (Gplus - Gminus) / (self.diff*vecl**2) - 2.0 * diag2 * vec / vecl**2
        grad = orthogopt(grad, self.coords, False)
        """
        C  Project out any component of the gradient along LOCALV (which is a unit vector)
        C  This is a big improvement for DFTB.
        """
        grad -= np.dot(grad, vec) * vec
        
        return diag2, grad
        
def testpot2():
    from pygmin.potentials.lj import LJ
    import itertools
    pot = LJ()
    a = 1.12 #2.**(1./6.)
    theta = 20./360*np.pi
    coords = [ 0., 0., 0., \
              -a, 0., 0., \
              a*np.cos(theta), a*np.sin(theta), 0. ]
    c = np.reshape(coords, [3,3])
    for i, j in itertools.combinations(range(3), 2):
        r = np.linalg.norm(c[i,:] - c[j,:])
        print i, j, r 

def testpot1():
    from pygmin.potentials.lj import LJ
    import itertools
    pot = LJ()
    a = 1.12 #2.**(1./6.)
    theta = 60./360*np.pi
    coords = [ 0., 0., 0., \
              -a, 0., 0., \
              -a/2, a*np.cos(theta), 0., \
              -a/2, -a*np.cos(theta), 0.1 \
              ]
    natoms = len(coords)/3
    c = np.reshape(coords, [-1,3])
    for i, j in itertools.combinations(range(natoms), 2):
        r = np.linalg.norm(c[i,:] - c[j,:])
        print i, j, r 
    
    e, g = pot.getEnergyGradient(coords)
    print "initial E", e
    print "initial G", g, np.linalg.norm(g)

    eigpot = LowestEigPot(coords, pot)
    vec = np.random.rand(len(coords))
    e, g = eigpot.getEnergyGradient(vec)
    print "eigenvalue", e 
    print "eigenvector", g
    
    if True:
        e, g, hess = pot.getEnergyGradientHessian(coords)
        print "shape hess", np.shape(hess)
        print "hessian", hess
        u, v = np.linalg.eig(hess)
        print "max imag value", np.max(np.abs(u.imag))
        print "max imag vector", np.max(np.abs(v.imag))
        u = u.real
        v = v.real
        print "eigenvalues", u
        for i in range(len(u)):
            print "eigenvalue", u[i], "eigenvector", v[:,i]
        #find minimum eigenvalue, vector
        imin = 0
        umin = 10.
        for i in range(len(u)):
            if np.abs(u[i]) < 1e-10: continue
            if u[i] < umin:
                umin = u[i]
                imin = i
        print "lowest eigenvalue ", umin, imin
        print "lowest eigenvector", v[:,imin]

    
    from pygmin.optimize.quench import lbfgs_py as quench
    ret = quench(vec, eigpot.getEnergyGradient, iprint=10, tol = 1e-5, maxstep = 1e-3, \
                 rel_energy = True)
    print ret
    
    print "lowest eigenvalue "
    print umin, imin
    print "lowest eigenvector"
    print v[:,imin]
    print "now the estimate"
    print ret[1]
    print ret[0]


if __name__ == "__main__":
    testpot1()
    
    
    