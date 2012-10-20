import numpy as np
import copy
from collections import namedtuple


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
    def __init__(self, pot, coords, verbose=False, orthogZeroEigs = 0):
        """
        :orthogZeroEigs: the function which makes a vector orthogonal to known zero 
            eigenvectors
            default value is 0, which means use the default function orthogopt.
            if None is pass then no function will be used
        """
        self.pot = pot
        self.eigvec = np.random.rand(len(coords)) #initial random guess
        self.verbose = verbose
        self.H0 = None
        self.orthogZeroEigs = orthogZeroEigs

    
    def getLowestEigenvalue(self, coords, iprint=400, tol=1e-6, nsteps=500, **kwargs):
        eigpot = LowestEigPot(coords, self.pot, orthogZeroEigs=self.orthogZeroEigs)
        from pygmin.optimize.lbfgs_py import LBFGS
        quencher = LBFGS(self.eigvec, eigpot, 
                         rel_energy=True, H0=self.H0, **kwargs)
        ret = quencher.run(iprint=iprint, tol=tol, nsteps=nsteps)
        self.H0 = quencher.H0
        #ret = quench(self.eigvec, eigpot.getEnergyGradient, iprint=400, tol = 1e-5, maxstep = 1e-3, rel_energy=True)
        self.eigval = ret[1]
        self.eigvec = ret[0]
        
        if False and self.verbose: #debugging
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

        #check if success
        rms = ret[2]
        return rms <= tol 

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



class FindTransitionState(object):
    def __init__(self, coords, pot, tol = 1e-4, event=None, nsteps=1000, 
                 tsSearchParams = None, nfail_max=5, **kwargs):
        """
        this is the class which implements the routine for finding the transition state
        """
        """
        notes: this class needs to deal with
        
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
            
        
        **tolerance for any of the minimizations must be at least as tight as the total minimization or it will never end
        
        note: This, when the lowest eigenvalue search is repeatedly failing, this routine can
            take a very long time.  We should recognize those  situations and fail early.
            Probably we should recognize when the lowest eigenvalue search fails and just end 
            if it fails 10 times in row.
        """
        self.pot = pot
        self.tangent_space_quencher = defaults.quenchRoutine
        self.coords = np.copy(coords)
        self.tol = tol
        self.nsteps = nsteps
        self.event = event
        self.nfail_max = nfail_max
        self.nfail = 0
        #self.tangent_space_quench_params = copy.copy(defaults.quenchParams)

        if tsSearchParams is None:
            tsSearchParams = defaults.tsSearchParams
        self.iprint = tsSearchParams.get("iprint")
        if self.iprint is None: iprint = -1
        if tsSearchParams.has_key("orthogZeroEigs"):
            has_orthogZeroEigs = True
            orthogZeroEigs = tsSearchParams["orthogZeroEigs"]  #this could meaningfully be None
        else:
            has_orthogZeroEigs = False
    
        self.lowestEigenvectorQuenchParams = defaults.lowestEigenvectorQuenchParams
            
                
        if has_orthogZeroEigs:
            self.tspot = TSRefinementPotential(pot, coords, orthogZeroEigs=orthogZeroEigs, **kwargs)
        else:
            self.tspot = TSRefinementPotential(pot, coords, **kwargs)
        
        self.rmsnorm = 1./np.sqrt(float(len(coords))/3.)
        self.oldeigvec = None

    def run(self):
        coords = np.copy(self.coords)
        for i in xrange(self.nsteps):
            
            #get the lowest eigenvalue and eigenvector
            coords = self.getLowestEigenVector(coords, i)
            
            #step uphill along the direction of the lowest eigenvector
            coords = self.stepUphill(coords)

            if False:
                #maybe we want to update the lowest eigenvector now that we've moved?
                #david thinks this is a bad idea
                coords = self.getLowestEigenVector(coords, i)

            #minimize the coordinates in the space perpendicular to the lowest eigenvector
            coords, rms = self.minimizeTangentSpace(coords)


            #check if we are done and print some stuff
            E, grad = self.pot.getEnergyGradient(coords)
            rms = np.linalg.norm(grad) * self.rmsnorm
            gradpar = np.dot(grad, self.tspot.eigvec) / np.linalg.norm(self.tspot.eigvec)
            
            if self.iprint > 0:
                if i % self.iprint == 0:
                    print "findTransitionState:", i, E, rms, "eigenvalue", self.tspot.eigval, "rms perpendicular", rms, "grad parallel", gradpar
            
            if callable(self.event):
                self.event(E, coords, rms)
            if rms < self.tol:
                break
            if self.nfail > self.nfail_max:
                print "findTransitionState fail"
                break

        #done, print some data
        print "findTransitionState done:", i, E, rms, "eigenvalue", self.tspot.eigval
    
        #check if results make sense
        if self.tspot.eigval >= 0.:
            print "warning: transition state has positive eigenvalue", self.tspot.eigval
        if rms > self.tol:
            print "warning: transition state search appears to have failed: rms", rms

        #return results
        return namedtuple("TransitionStateResults", "coords,energy,eigenval,eigenvec,grad,rms,nsteps")(
                    coords, E, self.tspot.eigval, self.tspot.eigvec, grad, rms, i)



        
    def getLowestEigenVector(self, coords, i):
        if self.lowestEigenvectorQuenchParams.has_key("tol"):
            tol = self.lowestEigenvectorQuenchParams["tol"]
        else:
            tol = 1e-6
            self.lowestEigenvectorQuenchParams["tol"] = tol

        success = self.tspot.getLowestEigenvalue(coords, **self.lowestEigenvectorQuenchParams)
        if self.tspot.eigval > 0.:
            print "warning transition state search found positive lowest eigenvalue", self.tspot.eigval, \
                "step", i
            if i == 0: 
                print "WARNING *** initial eigenvalue is positive - increase NEB spring constant?"
        if i > 0:
            overlap = np.dot(self.oldeigvec, self.tspot.eigvec)
            if overlap < 0.5:
                print "warning: the new eigenvector has low overlap with previous", overlap
        self.oldeigvec = self.tspot.eigvec.copy()

        if success:
            self.nfail = 0
        else:
            self.nfail += 1

        return coords
    
    def minimizeTangentSpace(self, coords):
        """
        now minimize the energy in the space perpendicular to eigvec.
        There's no point in spending much effort on this until 
        we've gotten close to the transition state.  So limit the number of steps
        to 10 until we get close.
        """
        tol = self.tol
        E, grad = self.pot.getEnergyGradient(coords)
        gradpar = np.dot(grad, self.tspot.eigvec) / np.linalg.norm(self.tspot.eigvec)
        nstepsperp = 10
        if np.abs(gradpar) <= tol*2.:
            nstepsperp = 100

        ret = self.tangent_space_quencher(coords, self.tspot.getEnergyGradient, nsteps=nstepsperp, tol=tol*0.2)
        coords = ret[0]
        return coords, ret[2]

    def stepUphill(self, coords):
        self.tspot.stepUphill(coords)
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
    neb.optimize(quenchParams={"iprint": 30, "nsteps":100})
    neb.MakeAllMaximaClimbing()
    neb.optimize(quenchParams={"iprint": 30, "nsteps":100})
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
        ret = findTransitionStateNew(coords, pot, event=printevent, verbose = False)
        
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
