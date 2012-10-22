import numpy as np
import copy
from collections import namedtuple
from scipy.optimize import Result

from lowest_eig_pot import LowestEigPot, findLowestEigenVector
from pygmin.potentials.potential import potential as basepot
from pygmin.storage.savenlowest import SaveN
import pygmin.defaults as defaults
from pygmin.defaults import lowestEigenvectorQuenchParams

__all__ = ["findTransitionState", "FindTransitionState"]

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


class TSRefinementPotential(basepot):
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
    def __init__(self, pot, eigenvec):
        """
        Parameters
        ----------
        
        :orthogZeroEigs: the function which makes a vector orthogonal to known zero 
            eigenvectors
            The default value is 0, which means use the default function orthogopt which assumes
            rotational and translational invariance.
            If None is pass then no function will be used
        """
        self.pot = pot
        self.eigenvec = eigenvec

    


    
    def getEnergyGradient(self, coords):
        """
        return the energy and the gradient with the component along the eigenvec removed.
        For use in energy minimization in the space perpendicular to eigenvec
        """
        e, grad = self.pot.getEnergyGradient(coords)
        #norm = np.sum(self.eigenvec)
        grad -= np.dot(grad, self.eigenvec) * self.eigenvec
        return e, grad



class FindTransitionState(object):
    """
    todo:
        if the eigenvalue sign goes from positive to negative,
        go back to where it was negative and take a smaller step
    """
    def __init__(self, coords, pot, tol=1e-4, event=None, nsteps=1000, 
                 nfail_max=5, eigenvec=None, iprint=-1, orthogZeroEigs=0,
                 lowestEigenvectorQuenchParams=dict()
                 ):
        """
        This class implements the routine for finding the nearest transition state
        
        Parameters
        ----------
        coords : 
            the starting coordinates
        pot : 
            the potential class
        tol : 
            the tolerance for the rms gradient
        event : callable
            This will be called after each step
        nsteps : 
            number of iterations
        nfail_max :
            if the lowest eigenvector search fails this many times in a row than the
            algorithm ends
        kwargs : 
            additional parameters passed to the
        eigenvec : 
            a guess for the initial lowest eigenvector
        
            
        
        Notes
        -----
        
        It is composed of the following steps
            1) Find eigenvector corresponding to the lowest *nonzero* eigenvector.  
            
            2) Step uphill in the direction of the lowest eigenvector
            
            3) minimize in the space tangent to the lowest eigenvector
         
        """
        """
        implementation notes: this class needs to deal with
        
        params for stepUphill : 
            probably only maxstep
        
        params for lowest eigenvector search : 
            should be passable and loaded from defaults.
            important params : 
                max steps 
                what to do if a negative eigenvector is found
                what to do if small overlap with previous eigenvector
        
        params for tangent space search : 
            should be passable and loaded from defaults.
            
        
        **tolerance for any of the minimizations must be at least as tight 
        as the total minimization or it will never end
        
        """
        self.pot = pot
        self.tangent_space_quencher = defaults.quenchRoutine
        self.coords = np.copy(coords)
        self.tol = tol
        self.nsteps = nsteps
        self.event = event
        self.nfail_max = nfail_max
        self.nfail = 0
        self.eigenvec = eigenvec
        self.orthogZeroEigs = orthogZeroEigs
        self.iprint = iprint
        self.lowestEigenvectorQuenchParams = lowestEigenvectorQuenchParams
        #self.tangent_space_quench_params = copy.copy(defaults.quenchParams)
            
        
        self.rmsnorm = 1./np.sqrt(float(len(coords))/3.)
        self.oldeigenvec = None

    def run(self):
        coords = np.copy(self.coords)
        for i in xrange(self.nsteps):
            
            #get the lowest eigenvalue and eigenvector
            self.getLowestEigenVector(coords, i)
            
            #step uphill along the direction of the lowest eigenvector
            coords = self.stepUphill(coords)

            if False:
                #maybe we want to update the lowest eigenvector now that we've moved?
                #david thinks this is a bad idea
                self.getLowestEigenVector(coords, i)

            #minimize the coordinates in the space perpendicular to the lowest eigenvector
            coords, tangentrms = self.minimizeTangentSpace(coords)


            #check if we are done and print some stuff
            E, grad = self.pot.getEnergyGradient(coords)
            rms = np.linalg.norm(grad) * self.rmsnorm
            gradpar = np.dot(grad, self.eigenvec) / np.linalg.norm(self.eigenvec)
            
            if self.iprint > 0:
                if i % self.iprint == 0:
                    print "findTransitionState:", i, E, rms, "eigenvalue", self.eigenval, "rms perpendicular", tangentrms, "grad parallel", gradpar
            
            if callable(self.event):
                self.event(E, coords, rms)
            if rms < self.tol:
                break
            if self.nfail >= self.nfail_max:
                print "stopping findTransitionState.  too many failures in eigenvector search"
                break

        #done.  do one last eigenvector search because coords may have changed
        self.getLowestEigenVector(coords, i)

        #done, print some data
        print "findTransitionState done:", i, E, rms, "eigenvalue", self.eigenval
    
        #check if results make sense
        if self.eigenval >= 0.:
            print "warning: transition state is ending with positive eigenvalue", self.eigenval
        if rms > self.tol:
            print "warning: transition state search appears to have failed: rms", rms
            success = False
        else:
            success = True

        #return results
        res = Result()
        res.coords = coords
        res.energy = E
        res.eigenval = self.eigenval
        res.eigenvec = self.eigenvec
        res.grad = grad
        res.rms = rms
        res.nsteps = i
        res.success = success
        return res


        
    def getLowestEigenVector(self, coords, i):
        if not hasattr(self, "H0"):
            self.H0 = None
        res = findLowestEigenVector(coords, self.pot, H0=self.H0, eigenvec0=self.eigenvec, 
                                    orthogZeroEigs=self.orthogZeroEigs,
                                    **self.lowestEigenvectorQuenchParams)
        
        if res.eigenval > 0.:
            print "warning transition state search found positive lowest eigenvalue", res.eigenval, \
                "step", i
            if i == 0: 
                print "WARNING *** initial eigenvalue is positive - increase NEB spring constant?"
        if i > 0:
            overlap = np.dot(self.oldeigenvec, res.eigenvec)
            if overlap < 0.5:
                print "warning: the new eigenvector has low overlap with previous", overlap
        
        self.H0 = res.H0
        self.eigenvec = res.eigenvec
        self.eigenval = res.eigenval
        self.oldeigenvec = self.eigenvec.copy()
        
        if res.success:
            self.nfail = 0
        else:
            self.nfail += 1

    
    def minimizeTangentSpace(self, coords):
        """
        now minimize the energy in the space perpendicular to eigenvec.
        There's no point in spending much effort on this until 
        we've gotten close to the transition state.  So limit the number of steps
        to 10 until we get close.
        """
        tol = self.tol
        E, grad = self.pot.getEnergyGradient(coords)
        gradpar = np.dot(grad, self.eigenvec) / np.linalg.norm(self.eigenvec)
        nstepsperp = 10
        if np.abs(gradpar) <= tol*2.:
            nstepsperp = 100

        tspot = TSRefinementPotential(self.pot, self.eigenvec)
        ret = self.tangent_space_quencher(coords, tspot.getEnergyGradient, nsteps=nstepsperp, tol=tol*0.2)
        coords = ret[0]
        rms = ret[2]
        return coords, rms

    def stepUphill(self, coords, maxstep = 0.1):
        """
        step uphill in the direction of self.eigenvec.  self.eigenval is used
        to determine the best stepsize
        """
        e, grad = self.pot.getEnergyGradient(coords)
        F = np.dot(grad, self.eigenvec) 
        h = 2.*F/ np.abs(self.eigenval) / (1. + np.sqrt(1.+4.*F**2/self.eigenval**2 ))

        if np.abs(h) > maxstep:
            h *= maxstep / abs(h)
        coords += h * self.eigenvec

        return coords

def findTransitionState(*args, **kwargs):
    """
    simply a wrapper for initializing and running FindTransitionState
    """
    finder = FindTransitionState(*args, **kwargs)
    return finder.run()



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
    #neb.optimize()
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
    neb.optimize(quenchParams={"iprint": -30, "nsteps":100})
    neb.MakeAllMaximaClimbing()
    #neb.optimize(quenchParams={"iprint": 30, "nsteps":100})
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

        defaults.quenchParams["iprint"] = 1
        #defaults.lowestEigenvectorQuenchParams["iprint"] = 1
        defaults.tsSearchParams["iprint"] = 1
        
        printevent = PrintEvent(fout)
        print ""
        print "starting the transition state search"
        ret = findTransitionState(coords, pot, event=printevent, iprint=-1)
        print ret
        #coords, eval, evec, e, grad, rms = ret
        e = pot.getEnergy(ret.coords)
        printxyz(fout, coords2, line2=str(e))

    print "finished searching for transition state"
    print "energy", e
    print "rms grad", ret.rms
    print "eigenvalue", ret.eigenval
    
    if False:
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
