import numpy as np
from lowest_eig_pot import LowestEigPot
from pygmin.potentials.potential import potential as basepot
from pygmin.storage.savenlowest import SaveN

def analyticalLowestEigenvalue(coords, pot):
    e, g, hess = pot.getEnergyGradientHessian(coords)
    #print "shape hess", np.shape(hess)
    #print "hessian", hess
    u, v = np.linalg.eig(hess)
    #print "max imag value", np.max(np.abs(u.imag))
    #print "max imag vector", np.max(np.abs(v.imag))
    u = u.real
    v = v.real
    #print "eigenvalues", u
    #for i in range(len(u)):
    #    print "eigenvalue", u[i], "eigenvector", v[:,i]
    #find minimum eigenvalue, vector
    imin = 0
    umin = 10.
    for i in range(len(u)):
        if np.abs(u[i]) < 1e-10: continue
        if u[i] < umin:
            umin = u[i]
            imin = i
    #print "analytical lowest eigenvalue ", umin, imin
    #print "lowest eigenvector", v[:,imin]
    return umin, v[:,imin]


class TransitionStateRefinement(basepot):
    """
    this is a potential wrapper for use in an optimization package.  The intent
    is to try to locate the nearest transition state of the system.
    
    For each call of getEnergyGradient it will
    
    1) try to estimate the lowest nonzero eigenvalue and corresponding eigenvector of the Hermitian
    
    2) return a gradient vector which points downhill in all directions, but uphill in the direction 
    of the eigenvector
    """
    def __init__(self, pot, coords):
        print "WARNING: TransitionStateRefinement is still in beta and probably not working"
        self.pot = pot
        self.eigvec = np.random.rand(len(coords))
        self.H0 = None    
    
    def getEnergyGradient(self, coords):
        
        #try to estimate the lowest nonzero eigenvalue and corresponding eigenvector of the Hermitian
        eigpot = LowestEigPot(coords, self.pot)
        from pygmin.optimize.quench import lbfgs_py as quench
        from pygmin.optimize.lbfgs_py import LBFGS
        print ""
        print "estimating lowest eigenvalue and eigenvector"
        quencher = LBFGS(self.eigvec, eigpot, maxstep=1e-3, \
                         rel_energy = True, H0 = self.H0)
        ret = quencher.run(iprint = 400, tol=1e-6)
        self.H0 = quencher.H0
        #ret = quench(self.eigvec, eigpot.getEnergyGradient, iprint=400, tol = 1e-5, maxstep = 1e-3, rel_energy=True)
        self.eigval = ret[1]
        self.eigvec = -ret[0]

        
        print "eigenvalue is estimated to be", self.eigval
        if True:
            trueval, truevec = analyticalLowestEigenvalue(coords, self.pot)
            print "analytical lowest eigenvalue", trueval
            maxdiff = np.max(np.abs(truevec - self.eigvec))
            print "maximum difference between estimated and analytical eigenvectors", maxdiff, \
                np.linalg.norm(self.eigvec), np.linalg.norm(truevec), np.dot(truevec, self.eigvec)
            if True:
                print self.eigvec
                print truevec
            
        
        
        self.e, self.grad = self.pot.getEnergyGradient(coords)
        self.rms = np.linalg.norm(self.grad) / np.sqrt(len(self.grad)/3)
        print "energy rms", self.e, self.rms
        #remove the component of the gradient along the eigenvector
        F = np.dot(self.grad, self.eigvec)
        
        #determine how large a step to take along the eigenvector
        h = 2.*F/ np.abs(self.eigval) / (1. + np.sqrt(1.+4.*F**2/self.eigval**2 ))
        print "stepsize h", h
        newgrad = self.grad + (- F + h) * self.eigvec
        
        print "rms of newgrad", np.linalg.norm(newgrad) / np.sqrt(len(newgrad)/3)
        print "dot(eigvec, grad)", np.dot(self.grad, self.eigvec)
        
        self.h = h
        
        return 0., newgrad

class NegativePot(basepot):
    def __init__(self, pot):
        self.pot = pot
    def getEnergy(self, coords):
        e = self.pot.getEnergy(coords)
        return -e

    def getEnergyGradient(self, coords):
        e, g = self.pot.getEnergyGradient(coords)
        return -e, -g

class TangentPot(basepot):
    """
    project the N dimensional potential into a one dimensional 
    problem along the line defined by X + a*V
    """
    def __init__(self, vec, pot):
        self.V = vec.copy()
        self.pot = pot
        self.Vnorm = np.linalg.norm(self.V)
        self.V /= self.Vnorm
    def getEnergy(self, coords):
        return self.pot.getEnergy( coords )
    def getEnegyGradient(self, coords):
        e, g = self.pot.getEnergyGradient( coords )
        #project g onto V
        g -= np.dot(self.V, g)*self.V
        return e, g

def findTransitionState(coords, pot):
    from pygmin.optimize.lbfgs_py import LBFGS
    from pygmin.optimize.bfgs import lineSearch as linesearch
    from pygmin.optimize.quench import lbfgs_py as quench
    tspot = TransitionStateRefinement(pot, coords)
    negpot = NegativePot(pot)

    
    lbfgs = LBFGS(coords, pot)
    neglbfgs = LBFGS(coords, negpot)
    
    nsteps = 100
    fout =open("out1.xyz", "w")
    from pygmin.printing.print_atoms_xyz import printAtomsXYZ as printxyz
    for i in xrange(nsteps):
        p, dx = tspot.getEnergyGradient(coords)
        E = tspot.e
        grad = tspot.grad
        eigvec = tspot.eigvec
        
                
        print "maximize the energy in the direction parallel to eigvec"
        if True:
            h = tspot.h
            if np.abs(h) > 0.1:
                h *= 0.1 / abs(h)
            dx = h * eigvec
            coords += h*eigvec
            e = pot.getEnergy(coords)
        else:
            coords, E, grad = neglbfgs.takeStepNoLineSearch(coords, -E, -grad, dx)
            E = -E
            grad = -grad

        
        print "minimize the energy in the direction perpendicular to eigvec"
        tangentpot = TangentPot(eigvec, pot)
        ret = quench(coords, tangentpot.getEnegyGradient)
        coords = ret[0]
        E = ret[1]
        
        printxyz(fout, coords, line2=str(E))


    fout.close()


def testgetcoordsLJ():
    a = 1.12 #2.**(1./6.)
    theta = 60./360*np.pi
    coords = [ 0., 0., 0., \
              -a, 0., 0., \
              -a/2, a*np.cos(theta), 0., \
              -a/2, -a*np.cos(theta), 0.3 \
              ]
    coords = np.array(coords)
    return coords


def guesstsATLJ():
    from pygmin.potentials.ATLJ import ATLJ
    pot = ATLJ(Z = 2.)
    a = 1.12 #2.**(1./6.)
    theta = 60./360*np.pi
    coords1 = np.array([ 0., 0., 0., \
              -a, 0., 0., \
              -a/2, -a*np.cos(theta), 0. ])
    coords2 = np.array([ 0., 0., 0., \
              -a, 0., 0., \
              a, 0., 0. ])
    from pygmin.optimize.quench import lbfgs_py as quench
    ret1 = quench(coords1, pot.getEnergyGradient)
    ret2 = quench(coords2, pot.getEnergyGradient)
    coords1 = ret1[0]
    coords2 = ret2[0]
    from pygmin.NEB.NEB import NEB
    neb = NEB(coords1, coords2, pot)
    neb.optimize()
    neb.MakeAllMaximaClimbing()
    neb.optimize()
    for i in xrange(len(neb.energies)):
        if(neb.isclimbing[i]):
            coords = neb.coords[i,:]
    return pot, coords

def guessts(coords1, coords2, pot):
    from pygmin.optimize.quench import lbfgs_py as quench
    ret1 = quench(coords1, pot.getEnergyGradient)
    ret2 = quench(coords2, pot.getEnergyGradient)
    coords1 = ret1[0]
    coords2 = ret2[0]
    from pygmin.NEB.NEB import NEB
    neb = NEB(coords1, coords2, pot)
    neb.optimize()
    neb.MakeAllMaximaClimbing()
    neb.optimize()
    for i in xrange(len(neb.energies)):
        if(neb.isclimbing[i]):
            coords = neb.coords[i,:]
    return pot, coords, neb.coords[0,:], neb.coords[-1,:]


def guesstsLJ():
    from pygmin.potentials.lj import LJ
    pot = LJ()
    natoms = 9
    coords = np.random.uniform(-1,1,natoms*3)
    from pygmin.basinhopping import BasinHopping
    from pygmin.takestep.displace import RandomDisplacement
    from pygmin.takestep.adaptive import AdaptiveStepsize
    from pygmin.storage.savenlowest import SaveN
    saveit = SaveN(10)
    takestep1 = RandomDisplacement()
    takestep = AdaptiveStepsize(takestep1, frequency=15)
    bh = BasinHopping(coords, pot, takestep, storage=saveit, outstream=None)
    bh.run(100)
    coords1 = saveit.data[0].coords
    coords2 = saveit.data[1].coords
    
    return guessts(coords1, coords2, pot)




    

def testgetcoordsATLJ():
    a = 1.12 #2.**(1./6.)
    theta = 40./360*np.pi
    coords = [ 0., 0., 0., \
              -a, 0., 0., \
              a*np.cos(theta), a*np.sin(theta), 0. ]
    return np.array(coords)

def testpot1():
    import itertools
    from pygmin.printing.print_atoms_xyz import printAtomsXYZ as printxyz
    pot, coords, coords1, coords2 = guesstsLJ()
    coordsinit = np.copy(coords)
    natoms = len(coords)/3
    c = np.reshape(coords, [-1,3])
    for i, j in itertools.combinations(range(natoms), 2):
        r = np.linalg.norm(c[i,:] - c[j,:])
        print i, j, r 
    
    e, g = pot.getEnergyGradient(coords)
    print "initial E", e
    print "initial G", g, np.linalg.norm(g)
    

    tspot = TransitionStateRefinement(pot, coords)
    #e, g = tspot.getEnergyGradient(coords)
    print "so called energy", e 
    print "gradient", g
    
    
    from pygmin.printing.print_atoms_xyz import PrintEvent
    from pygmin.optimize.quench import steepest_descent as quench
    #ret = quench(coords, tspot.getEnergyGradient, iprint=1, tol = 1e-2, maxstep = 1e-1, dx = 1e-2)
    #ret = quench(coords, tspot.getEnergyGradient, iprint=1, tol = 1e-2, maxstep = 1e-1)

    #print ret
    
    with open("out.xyz", "w") as fout:
        e = pot.getEnergy(coords1)
        printxyz(fout, coords1, line2=str(e))
        e = pot.getEnergy(coordsinit)
        printxyz(fout, coordsinit, line2=str(e))
        e = pot.getEnergy(coords2)
        printxyz(fout, coords2, line2=str(e))
        
        findTransitionState(coords, pot)

        if False:
            printevent = PrintEvent(fout)
            ret = quench(coords, tspot.getEnergyGradient, iprint=1, tol = 1e-3, maxstep = 1e-3, \
                        dx = 1e-4, nsteps=100, event=printevent)
            e = pot.getEnergy(ret[0])
            printxyz(fout, ret[0], line2=str(e))

    


if __name__ == "__main__":
    testpot1()
