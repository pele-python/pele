import numpy as np
import copy
import logging

from pele.optimize import Result, MYLBFGS
from pele.optimize import mylbfgs
from pele.potentials.potential import BasePotential
from pele.transition_states import findLowestEigenVector


__all__ = ["findTransitionState", "FindTransitionState"]

logger = logging.getLogger("pele.connect.findTS")

class _TransversePotential(BasePotential):
    """This wraps a potential and returns the gradient with the component parallel to a vector removed
    
    Parameters
    ----------
    
    
        
    """
    def __init__(self, potential, vector):
        """
        """
        self.pot = potential
        self.update_vector(vector)
        self.nfev = 0
    
    def update_vector(self, vector):
        self.vector = vector.copy()
        # normalize
        self.vector /= np.linalg.norm(vector)

    def projected_energy_gradient(self, true_energy, true_gradient):
        grad = true_gradient - np.dot(true_gradient, self.vector) * self.vector
        return true_energy, grad

        
    def getEnergyGradient(self, coords):
        """
        return the energy and the gradient with the component along the
        eigenvec removed.  For use in energy minimization in the space
        perpendicular to eigenvec
        """
        self.nfev += 1
        e, grad = self.pot.getEnergyGradient(coords)
        self.true_gradient = grad.copy()
        self.true_energy = e
        #norm = np.sum(self.eigenvec)
        return self.projected_energy_gradient(e, grad)



class FindTransitionState(object):
    """
    This class implements the hybrid eigenvector following routine for finding the nearest transition state
    
    ***orthogZeroEigs is system dependent, don't forget to set it***
    
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
        if the lowest eigenvector search fails this many times in a row
        than the algorithm ends
    eigenvec0 : 
        a guess for the initial lowest eigenvector
    iprint :
        the interval at which to print status messages
    orthogZeroEigs : callable
        this function makes a vector orthogonal to the known zero
        eigenvectors

            orthogZeroEigs=0  : default behavior, assume translational and
                                rotational symmetry
            orthogZeroEigs=None : the vector is unchanged

    lowestEigenvectorQuenchParams : dict 
        these parameters are passed to the quench routine for he lowest
        eigenvector search 
    tangentSpaceQuenchParams : dict 
        these parameters are passed quench routine for the minimization in
        the space tabgent to the lowest eigenvector 
    max_uphill_step : 
        the maximum step uphill along the direction of the lowest
        eigenvector
    demand_initial_negative_vec : bool
        if True, abort if the initial lowest eigenvalue is positive
    negatives_before_check : int
        if start with positive eigenvector and demand_initial_negative_vec is False,
        the check to make sure that the eigenvalue is enabled after having had so 
        many negative eigenvalues before
    nsteps_tangent1, nsteps_tangent2 : int
        the number of iterations for tangent space minimization before and after
        the eigenvalue is deemed to be converged
    verbosity : int
        how much debugging information to print (only partially implemented)
        
    
    Notes
    -----
    
    It is composed of the following steps
        1) Find eigenvector corresponding to the lowest *nonzero*
        eigenvector.  
        
        2) Step uphill in the direction of the lowest eigenvector
        
        3) minimize in the space tangent to the lowest eigenvector
        
    The tolerances for the various steps of this algorithm must be correlated.
    if the tol for tangent space search is lower than the total tol, then it will never finish
    
    See Also
    --------
    findTransitionState : function wrapper for this class
    findLowestEigenVector : a core algorithm
    pele.landscape.LocalConnect : the class which most often calls this routine
    """
    def __init__(self, coords, pot, tol=1e-4, event=None, nsteps=100, 
                 nfail_max=200, eigenvec0=None, iprint=-1, orthogZeroEigs=0,
                 lowestEigenvectorQuenchParams=dict(),
                 tangentSpaceQuenchParams=dict(), 
                 max_uphill_step=0.1,
                 demand_initial_negative_vec=True,
                 negatives_before_check = 10,
                 nsteps_tangent1=10,
                 nsteps_tangent2=100,
                 verbosity=1,
                 first_order=False,
                 ):
        self.pot = pot
        self.coords = np.copy(coords)
        self.tol = tol
        self.nsteps = nsteps
        self.event = event
        self.nfail_max = nfail_max
        self.nfail = 0
        self.eigenvec = eigenvec0
        self.orthogZeroEigs = orthogZeroEigs
        self.iprint = iprint
        self.lowestEigenvectorQuenchParams = lowestEigenvectorQuenchParams
        self.max_uphill_step = max_uphill_step
        self.verbosity = verbosity
        self.tangent_space_quencher = mylbfgs #  should make this passable
        self.tangent_space_quench_params = dict(tangentSpaceQuenchParams.items())
        self.demand_initial_negative_vec = demand_initial_negative_vec    
        self.npositive_max = max(10, self.nsteps / 5)
        self.first_order = first_order
        
        self.rmsnorm = 1./np.sqrt(float(len(coords)))
        self.oldeigenvec = None

        #set tolerance for the tangent space minimization.  
        #Be sure it is tighter tolerance than self.tol
        self.tol_tangent = self.tol * 0.2
        if self.tangent_space_quench_params.has_key("tol"):
            self.tol_tangent = min(self.tol_tangent, 
                                   self.tangent_space_quench_params["tol"])
            del self.tangent_space_quench_params["tol"]
        self.nsteps_tangent1 = nsteps_tangent1
        self.nsteps_tangent2 = nsteps_tangent2
        if self.tangent_space_quench_params.has_key("maxstep"):
            self.maxstep_tangent = self.tangent_space_quench_params["maxstep"]
            del self.tangent_space_quench_params["maxstep"]
        else:            
            self.maxstep_tangent = 0.1 #this should be determined in a better way
        
        if not self.tangent_space_quench_params.has_key("logger"):
            self.tangent_space_quench_params["logger"] = logging.getLogger("pele.connect.findTS.tangent_space_quench")


        #set some parameters used in finding lowest eigenvector
        #initial guess for Hermitian
        self.H0_leig = None 
        
        self.H0_transverse = None
        
        self.reduce_step = 0
        self.step_factor = .1
        self.npositive = 0
        
    @classmethod
    def params(cls, obj = None):
        if obj is None:
            obj = FindTransitionState(np.zeros(2), None)
        
        params = dict()
        
        params["tangentSpaceQuenchParams"] = obj.tangent_space_quench_params.copy()
        params["lowestEigenvectorQuenchParams"] = obj.lowestEigenvectorQuenchParams.copy()
        params["tol"] = obj.tol
        params["nsteps"] = obj.nsteps
        params["nfail_max"] = obj.nfail_max
        params["iprint"] = obj.iprint
        params["max_uphill_step"]=obj.max_uphill_step
        params["demand_initial_negative_vec"]=obj.demand_initial_negative_vec
        params["nsteps_tangent1"]=obj.nsteps_tangent1
        params["nsteps_tangent2"]=obj.nsteps_tangent2
        params["verbosity"]=obj.verbosity
        
        # event=None, eigenvec0=None, orthogZeroEigs=0,
        return params
    
    
    
    def _saveState(self, coords):
        self.saved_coords = np.copy(coords)
        self.saved_eigenvec = np.copy(self.eigenvec)
        self.saved_eigenval = self.eigenval
        self.saved_overlap = self.overlap
        self.saved_H0_leig = self.H0_leig
        self.saved_H0_transverse = self.H0_transverse
        self.saved_energy = self.energy
        self.saved_gradient = self.gradient.copy()
        #self.saved_oldeigenvec = np.copy(self.oldeigenvec)

    def _resetState(self):
        coords = np.copy(self.saved_coords)
        self.eigenvec = np.copy(self.saved_eigenvec)
        self.eigenval = self.saved_eigenval
        self.oldeigenvec = np.copy(self.eigenvec)
        self.overlap = self.saved_overlap
        self.H0_leig = self.saved_H0_leig
        self.H0_transverse = self.saved_H0_transverse
        self.energy = self.saved_energy
        self.gradient = self.saved_gradient.copy()
        return coords

    def _compute_gradients(self, coords):
#        self._transverse_energy, self._transverse_gradient = self.transverse_potential.getEnergyGradient(coords)
#        self.energy = self.transverse_potential.true_energy
#        self.gradient = self.transverse_potential.true_gradient.copy()
        self.energy, self.gradient = self.pot.getEnergyGradient(coords)

    def _set_energy_gradient(self, energy, gradient):
        self.energy = energy
        self.gradient = gradient.copy()

    def get_energy(self):
        return self.energy
    
    def get_gradient(self):
        return self.gradient

    def run(self):
        """The main loop of the algorithm"""
        coords = np.copy(self.coords)
        res = Result() #  return object
        res.message = []
        
        # if starting with positive curvature, disable negative eigenvalue check
        # this will be reenabled as soon as the eigenvector becomes negative
        negative_before_check =  10

        self._compute_gradients(coords)
        for i in xrange(self.nsteps):
            
            # get the lowest eigenvalue and eigenvector
            self.overlap = self._getLowestEigenVector(coords, i)
            overlap = self.overlap
            #print self.eigenval
            
            if self.eigenval < 0:
                negative_before_check -= 1
            
            # check to make sure the eigenvector is ok
            if i == 0 or self.eigenval <= 0 or (negative_before_check > 0 and not self.demand_initial_negative_vec):
                self._saveState(coords)
                self.reduce_step = 0
            else:
                self.npositive += 1
                if self.npositive > self.npositive_max:
                    logger.warning( "positive eigenvalue found too many times. ending %s", self.npositive)
                    res.message.append( "positive eigenvalue found too many times %d" % self.npositive )
                    break
                if self.verbosity > 2:
                    logger.info("the eigenvalue turned positive. %s %s", self.eigenval, "Resetting last good values and taking smaller steps")
                coords = self._resetState()
                self.reduce_step += 1
            
            # step uphill along the direction of the lowest eigenvector
            coords = self._stepUphill(coords)
            self._compute_gradients(coords)

            # minimize the coordinates in the space perpendicular to the lowest eigenvector
            tangent_ret = self._minimizeTangentSpace(coords, energy=self.get_energy(), gradient=self.get_gradient())
            coords = tangent_ret.coords
            tangentrms = tangent_ret.rms


            # check if we are done and print some stuff
#            self._compute_gradients(coords) # this is unnecessary
            E = self.get_energy()
            grad = self.get_gradient()
            rms = np.linalg.norm(grad) * self.rmsnorm
            gradpar = np.dot(grad, self.eigenvec) / np.linalg.norm(self.eigenvec)
            
            if self.iprint > 0:
                if (i+1) % self.iprint == 0:
                    ostring = "findTS: %3d E %9g rms %8g eigenvalue %9g rms perp %8g grad par %9g overlap %g" % (
                                    i, E, rms, self.eigenval, tangentrms, gradpar, overlap)
                    extra = "  Evec search: %d rms %g" % (self.leig_result.nfev, self.leig_result.rms)
                    extra += "  Tverse search: %d step %g" % (self.tangent_result.nfev, 
                                                                    self.tangent_move_step)
                    extra += "  Uphill step:%g" % (self.uphill_step_size,)
                    logger.info("%s %s", ostring, extra)
            
            if callable(self.event):
                self.event(energy=E, coords=coords, rms=rms, eigenval=self.eigenval, stepnum=i)
            if rms < self.tol:
                break
            if self.nfail >= self.nfail_max:
                logger.warning("stopping findTransitionState.  too many failures in eigenvector search %s", self.nfail)
                res.message.append( "too many failures in eigenvector search %d" % self.nfail )
                break

            if i == 0 and self.eigenval > 0.:
                logger.warning("initial eigenvalue is positive - increase NEB spring constant?")
                if self.demand_initial_negative_vec:
                    logger.warning("            aborting transition state search")
                    res.message.append( "initial eigenvalue is positive %f" % self.eigenval )
                    break

        # done.  do one last eigenvector search because coords may have changed
        self._getLowestEigenVector(coords, i)

        # print some data
        logger.info("findTransitionState done: %s %s %s %s %s", i, E, rms, "eigenvalue", self.eigenval)
    
        success = True
        # check if results make sense
        if self.eigenval >= 0.:
            if self.verbosity > 2:
                logger.info( "warning: transition state is ending with positive eigenvalue %s", self.eigenval)
            success = False
        if rms > self.tol:
            if self.verbosity > 2:
                logger.info("warning: transition state search appears to have failed: rms %s", rms)
            success = False
        if i >= self.nsteps:
            res.message.append( "maximum iterations reached %d" % i )
            

        #return results
        res.coords = coords
        res.energy = E
        res.eigenval = self.eigenval
        res.eigenvec = self.eigenvec
        res.grad = grad
        res.rms = rms
        res.nsteps = i
        res.success = success
        return res


        
    def _getLowestEigenVector(self, coords, i, gradient=None):
        res = findLowestEigenVector(coords, self.pot, H0=self.H0_leig, eigenvec0=self.eigenvec, 
                                    orthogZeroEigs=self.orthogZeroEigs, first_order=self.first_order,
                                    gradient=None,
                                    **self.lowestEigenvectorQuenchParams)
        self.leig_result = res
        
#        if res.eigenval > 0.:
#            print "warning transition state search found positive lowest eigenvalue", res.eigenval, \
#                "step", i
        
        self.H0_leig = res.H0
        self.eigenvec = res.eigenvec
        self.eigenval = res.eigenval
        
        if i > 0:
            overlap = np.dot(self.oldeigenvec, res.eigenvec)
            if overlap < 0.5 and self.verbosity > 2:
                logger.info("warning: the new eigenvector has low overlap with previous %s %s", overlap, self.eigenval)
        else:
            overlap = 0.

        if res.success:
            self.nfail = 0
        else:
            self.nfail += 1
        
        self.oldeigenvec = self.eigenvec.copy()
        return overlap
    
    def _minimizeTangentSpace(self, coords, energy=None, gradient=None):
        """
        now minimize the energy in the space perpendicular to eigenvec.
        There's no point in spending much effort on this until 
        we've gotten close to the transition state.  So limit the number of steps
        to 10 until we get close.
        """
        #determine the number of steps
        #i.e. if the eigenvector is deemed to have converged
        use_gradpar = False
        if use_gradpar:
            E, grad = self.pot.getEnergyGradient(coords)
            gradpar = np.dot(grad, self.eigenvec) / np.linalg.norm(self.eigenvec)
            eigenvec_converged = np.abs(gradpar) <= self.tol*2.
        else:
            eigenvec_converged = self.overlap > .999 
        
        nstepsperp = self.nsteps_tangent1
        if eigenvec_converged:
            nstepsperp = self.nsteps_tangent2

        maxstep = self.maxstep_tangent
        if self.reduce_step > 0:
            maxstep *= (self.step_factor)**self.reduce_step



        tspot = _TransversePotential(self.pot, self.eigenvec)
        transverse_energy, transverse_gradient = tspot.projected_energy_gradient(energy, gradient) 
        coords1 = np.copy(coords)
        ret = self.tangent_space_quencher(coords, tspot, 
                                          nsteps=nstepsperp, tol=self.tol_tangent,
                                          maxstep=maxstep,
                                          H0 = self.H0_transverse,
                                          energy=transverse_energy, gradient=transverse_gradient,
                                          **self.tangent_space_quench_params)
        coords = ret.coords
        self.tangent_move_step = np.linalg.norm(coords - coords1)
        rms = ret.rms
        self.tangent_result = ret
        self.H0_transverse = self.tangent_result.H0
        self.energy = ret.energy
        try:
            self.gradient = tspot.true_gradient
        except AttributeError:
            # tspot was never called, use the same gradient
            if gradient is None:
                self._compute_gradients(coords)
            else:
                self.gradient = gradient
        return ret

    def _stepUphill(self, coords):
        """
        step uphill in the direction of self.eigenvec.  self.eigenval is used
        to determine the best stepsize
        """
        # the energy and gradient are already known
        e = self.get_energy()
        grad = self.get_gradient()
        F = np.dot(grad, self.eigenvec) 
        h = 2.*F/ np.abs(self.eigenval) / (1. + np.sqrt(1.+4.*F**2/self.eigenval**2 ))

        maxstep = self.max_uphill_step
        if self.reduce_step > 0:
            maxstep *= (self.step_factor)**self.reduce_step


        if np.abs(h) > maxstep:
            if self.verbosity >= 5:
                logger.debug("reducing step from %s %s %s", h, "to", maxstep) 
            h *= maxstep / abs(h)
        self.uphill_step_size = h
        coords += h * self.eigenvec

        return coords


def findTransitionState(*args, **kwargs):
    """
    simply a wrapper for initializing and running FindTransitionState
    
    See Also
    --------
    FindTransitionState : for all documentation
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
    from pele.potentials.ATLJ import ATLJ
    pot = ATLJ(Z = 2.)
    a = 1.12 #2.**(1./6.)
    theta = 60./360*np.pi
    coords1 = np.array([ 0., 0., 0., \
              -a, 0., 0., \
              -a/2, -a*np.cos(theta), 0. ])
    coords2 = np.array([ 0., 0., 0., \
              -a, 0., 0., \
              a, 0., 0. ])
    from pele.optimize import lbfgs_py as quench
    from pele.transition_states import InterpolatedPath
    ret1 = quench(coords1, pot.getEnergyGradient)
    ret2 = quench(coords2, pot.getEnergyGradient)
    coords1 = ret1[0]
    coords2 = ret2[0]
    from pele.transition_states import NEB
    neb = NEB(InterpolatedPath(coords1, coords2, 30), pot)
    neb.optimize()
    neb.MakeAllMaximaClimbing()
    #neb.optimize()
    for i in xrange(len(neb.energies)):
        if(neb.isclimbing[i]):
            coords = neb.coords[i,:]
    return pot, coords

def guessts(coords1, coords2, pot):
    from pele.optimize import lbfgs_py as quench
#    from pele.mindist.minpermdist_stochastic import minPermDistStochastic as mindist
    from pele.transition_states import NEB
    from pele.systems import LJCluster
    ret1 = quench(coords1, pot.getEnergyGradient)
    ret2 = quench(coords2, pot.getEnergyGradient)
    coords1 = ret1[0]
    coords2 = ret2[0]
    natoms = len(coords1)/3
    system = LJCluster(natoms)
    mindist = system.get_mindist()
    dist, coords1, coords2 = mindist(coords1, coords2) 
    print "dist", dist
    print "energy coords1", pot.getEnergy(coords1)
    print "energy coords2", pot.getEnergy(coords2)
    from pele.transition_states import InterpolatedPath
    neb = NEB(InterpolatedPath(coords1, coords2, 20), pot)
    #neb.optimize(quenchParams={"iprint" : 1})
    neb.optimize(iprint=-30, nsteps=100)
    neb.MakeAllMaximaClimbing()
    #neb.optimize(quenchParams={"iprint": 30, "nsteps":100})
    for i in xrange(len(neb.energies)):
        if(neb.isclimbing[i]):
            coords = neb.coords[i,:]
    return pot, coords, neb.coords[0,:], neb.coords[-1,:]


def guesstsLJ():
    from pele.potentials.lj import LJ
    pot = LJ()
    natoms = 9
    coords = np.random.uniform(-1,1,natoms*3)
    from pele.basinhopping import BasinHopping
    from pele.takestep.displace import RandomDisplacement
    from pele.takestep.adaptive import AdaptiveStepsize
    from pele.storage.savenlowest import SaveN
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
    from pele.printing.print_atoms_xyz import printAtomsXYZ as printxyz
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
    

    
    
    from pele.printing.print_atoms_xyz import PrintEvent

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
        ret = findTransitionState(coords, pot, iprint=-1)
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
        from pele.NEB.dimer import findTransitionState as dimerfindTS
        coords = coordsinit.copy()
        tau = np.random.uniform(-1,1,len(coords))
        tau /= np.linalg.norm(tau)
        ret = dimerfindTS(coords, pot, tau )
        enew, grad = pot.getEnergyGradient(ret.coords)
        print "energy", enew
        print "rms grad", np.linalg.norm(grad) / np.sqrt(float(len(ret.coords))/3.)



if __name__ == "__main__":
    testpot1()
