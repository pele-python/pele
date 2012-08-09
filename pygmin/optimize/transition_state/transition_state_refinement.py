import numpy as np
from lowest_eig_pot import LowestEigPot
from orthogopt import orthogopt
from pygmin.potentials.potential import potential as basepot
from pygmin.storage.savenlowest import SaveN
import pygmin.defaults as defaults

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
    
    Usage:
    
    getLowestEigenvalue():
        will attempt to calculate the lowest eigenvalue and corresponding eigenvector.
        
    stepUphill():
        will take a step uphill along the eigenvector
    
    getEnergyGradient():
        will return the gradient with the component along the eigenvector removed.  This is for
        energy minimization in the space tangent to the gradient    
    """
    def __init__(self, pot, coords, verbose=False, orthogZeroEigs = orthogopt):
        self.pot = pot
        self.eigvec = np.random.rand(len(coords)) #initial random guess
        self.verbose = verbose
        self.H0 = None
        self.orthogZeroEigs = orthogZeroEigs
    
    def getLowestEigenvalue(
                self, coords, iprint=400, tol = 1e-6, **kwargs):
        eigpot = LowestEigPot(coords, self.pot, orthogZeroEigs = self.orthogZeroEigs)
        from pygmin.optimize.lbfgs_py import LBFGS
        quencher = LBFGS(self.eigvec, eigpot, maxstep=1e-2, \
                         rel_energy = True, H0 = self.H0, **kwargs)
        ret = quencher.run(iprint = iprint, tol=tol)
        self.H0 = quencher.H0
        #ret = quench(self.eigvec, eigpot.getEnergyGradient, iprint=400, tol = 1e-5, maxstep = 1e-3, rel_energy=True)
        self.eigval = ret[1]
        self.eigvec = ret[0]
        
        if self.verbose: #debugging
            print ""
            print "eigenvalue is estimated to be", self.eigval
            try:
                trueval, truevec = analyticalLowestEigenvalue(coords, self.pot)
                print "analytical lowest eigenvalue", trueval
                print "overlap between estimated and analytical eigenvectors", \
                    np.dot(truevec, self.eigvec)
                if True:
                    print self.eigvec
                    print truevec
            except:
                pass

        
        return self.eigval, self.eigvec

    def stepUphill(self, coords, maxstep = 0.1):
        e, grad = self.pot.getEnergyGradient(coords)
        F = np.dot(grad, self.eigvec) 
        h = 2.*F/ np.abs(self.eigval) / (1. + np.sqrt(1.+4.*F**2/self.eigval**2 ))

        if np.abs(h) > maxstep:
            h *= maxstep / abs(h)
        coords += h * self.eigvec

    
    def getEnergyGradient(self, coords):
        """
        return the energy and the gradient with the component along the eigvec removed.
        For use in energy minimization in the space perpendicular to eigvec
        """
        e, grad = self.pot.getEnergyGradient(coords)
        #norm = np.sum(self.eigvec)
        grad -= np.dot(grad, self.eigvec) * self.eigvec
        return e, grad


def findTransitionState(coords, pot, tol = 1e-4, event=None, nsteps=1000, tsSearchParams = None, **kwargs):
    #from pygmin.optimize.quench import lbfgs_py as quench
    quenchRoutine = defaults.quenchRoutine
    tspot = TransitionStateRefinement(pot, coords, **kwargs)
    rmsnorm = 1./np.sqrt(float(len(coords))/3.)
    oldeigvec = None
    
    iprint = defaults.quenchParams.get("iprint")
    if iprint is None: iprint = -1

    for i in xrange(nsteps):
        tspot.getLowestEigenvalue(coords)
        if tspot.eigval > 0.:
            print "warning transition state search found positive lowest eigenvalue", tspot.eigval, \
                "step", i
        if i > 0:
            overlap = np.dot(oldeigvec, tspot.eigvec)
            if overlap < 0.5:
                print "warning: the new eigenvector has low overlap with previous", overlap
        oldeigvec = tspot.eigvec.copy()  
                
        #print "step uphill in the energy in the direction parallel to eigvec"
        tspot.stepUphill(coords)
        
        if True:
            #refine the eigenector again
            tspot.getLowestEigenvalue(coords)
            if tspot.eigval > 0.:
                print "warning transition state search found positive lowest eigenvalue", tspot.eigval, \
                    "step", i
            if i > 0:
                overlap = np.dot(oldeigvec, tspot.eigvec)
                if overlap < 0.5:
                    print "warning: the new eigenvector has low overlap with previous", overlap
            oldeigvec = tspot.eigvec.copy()  

        
        #get the component of the gradient parallel to eigvec
        E, grad = pot.getEnergyGradient(coords)
        gradpar = np.dot(grad, tspot.eigvec) / np.linalg.norm(tspot.eigvec)
        
        
        #print "minimize the energy in the direction perpendicular to eigvec"
        """
        now minimize the energy in the space perpendicular to eigvec.
        There's no point in spending much effort on this until 
        we've gotten close to the transition state.  So limit the number of steps
        to 10 until we get close.
        """
        nstepsperp = 10
        if np.abs(gradpar) <= tol*2.:
            nstepsperp = 1000
        ret = quenchRoutine(coords, tspot.getEnergyGradient, nsteps=nstepsperp, tol=tol*0.2)
        coords = ret[0]
        E, grad = pot.getEnergyGradient(coords)
        rms = np.linalg.norm(grad) * rmsnorm
        
        if iprint > 0:
            if i % iprint == 0:
                print "findTransitionState:", i, E, rms, "eigenvalue", tspot.eigval, "rms perpendicular", ret[2], "grad parallel", gradpar
        
        if event != None:
            event(E, coords, rms)
        if rms < tol:
            break
    
    if tspot.eigval >= 0.:
        print "warning: transition state has positive eigenvalue", tspot.eigval
    if rms > tol:
        print "warning: transition state search appears to have failed: rms", rms
    
    from collections import namedtuple
    return namedtuple("TransitionStateResults", "coords,energy,eigenval,eigenvec,grad,rms")(coords, E, tspot.eigval, tspot.eigvec, grad, rms)

###################################################################
#below here only stuff for testing
###################################################################

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
    from pygmin.mindist.minpermdist_stochastic import minPermDistStochastic as mindist
    from pygmin.NEB.NEB import NEB
    ret1 = quench(coords1, pot.getEnergyGradient)
    ret2 = quench(coords2, pot.getEnergyGradient)
    coords1 = ret1[0]
    coords2 = ret2[0]
    natoms = len(coords1)/3
    dist, coords1, coords2 = mindist(coords1, coords2, permlist=[range(natoms)])
    print "dist", dist
    print "energy coords1", pot.getEnergy(coords1)
    print "energy coords2", pot.getEnergy(coords2)
    neb = NEB(coords1, coords2, pot)
    #neb.optimize(quenchParams={"iprint" : 1})
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
    print ""
    

    
    
    from pygmin.printing.print_atoms_xyz import PrintEvent

    #print ret
    
    with open("out.xyz", "w") as fout:
        e = pot.getEnergy(coords1)
        print "energy of minima 1", e
        printxyz(fout, coords1, line2=str(e))
        e, grad = pot.getEnergyGradient(coordsinit)
        print "energy of NEB guess for the transition state", e, "rms grad", \
            np.linalg.norm(grad) / np.sqrt(float(len(coords))/3.)
        printxyz(fout, coordsinit, line2=str(e))
        e = pot.getEnergy(coords2)
        print "energy of minima 2", e
        printxyz(fout, coords2, line2=str(e))
        
        #mess up coords a bit
        coords += np.random.uniform(-1,1,len(coords))*0.05
        e = pot.getEnergy(coords)
        printxyz(fout, coords, line2=str(e))

        
        printevent = PrintEvent(fout)
        print ""
        print "starting the transition state search"
        ret = findTransitionState(coords, pot, event=printevent, verbose = False)
        
        #coords, eval, evec, e, grad, rms = ret
        e = pot.getEnergy(ret.coords)
        printxyz(fout, coords2, line2=str(e))

    print "finished searching for transition state"
    print "energy", e
    print "rms grad", ret.rms
    print "eigenvalue", ret.eigenval
    
    if True:
        print "now try the same search with the dimer method"
        from pygmin.NEB.dimer import findTransitionState as dimerfindTS
        coords = coordsinit.copy()
        tau = np.random.uniform(-1,1,len(coords))
        tau /= np.linalg.norm(tau)
        ret = dimerfindTS(coords, pot, tau )
        enew, grad = pot.getEnergyGradient(ret.coords)
        print "energy", enew
        print "rms grad", np.linalg.norm(grad) / np.sqrt(float(len(ret.coords))/3.)



if __name__ == "__main__":
    testpot1()
