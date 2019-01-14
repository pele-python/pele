from __future__ import print_function
import numpy as np
import logging

from pele.optimize import Result
from pele.optimize import mylbfgs
from pele.transition_states import FindLowestEigenVector
from pele.transition_states._dimer_translator import _DimerTranslator
from pele.transition_states._transverse_walker import _TransverseWalker
from pele.utils.hessian import get_smallest_eig

__all__ = ["findTransitionState", "FindTransitionState"]

logger = logging.getLogger("pele.connect.findTS")


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
    nsteps : 
        number of iterations
    eigenvec0 : 
        a guess for the initial lowest eigenvector
    nsteps_tangent1, nsteps_tangent2 : int
        the number of iterations for tangent space minimization before and after
        the eigenvalue is deemed to be converged
    event : callable
        This will be called after each step
    iprint :
        the interval at which to print status messages
    verbosity : int
        how much debugging information to print (only partially implemented)
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
    max_uphill_step_initial : float
        The initial maximum uphill step along the direction of the lowest 
        eigenvector.  The maximum uphill step is adjusted dynamically using a 
        trust radius
    max_uphill_step : float 
        The maximum value of the maximum uphill step. The maximum uphill step 
        is adjusted using a trust radius, but it can never become larger than 
        this.  
    check_negative : bool
        If True then the sign of the lowest eigenvalue is required to remain negative.
        If the sign becomes positive, then the step is retaken with smaller step size
    demand_initial_negative_vec : bool
        if True then abort if the initial lowest eigenvalue is positive
    negatives_before_check : int
        If starting with positive curvature, disable negative eigenvalue check.
        This will be re-enabled as soon as the eigenvalue becomes negative.
        Allow a certain number of good iterations with negative eigenvalues before
        demanding that the eigenvalue stay negative.
    nfail_max :
        if the lowest eigenvector search fails this many times in a row
        than the algorithm ends
    hessian_diagonalization : bool
        Diagonalize the Hessian matrix to find the lowest eigenvector rather
        than using the iterative procedure
        
    
    Notes
    -----
    
    It is composed of the following steps:
    
        1) Find eigenvector corresponding to the lowest *nonzero*
           eigenvector.  
        
        2) Step uphill in the direction of the lowest eigenvector
        
        3) minimize in the space tangent to the lowest eigenvector
        
    The tolerances for the various steps of this algorithm must be correlated.
    if the tolerance for tangent space search is lower than the total tolerance, 
    then it will never finish
    
    See Also
    --------
    findTransitionState : function wrapper for this class
    findLowestEigenVector : a core algorithm
    pele.landscape.LocalConnect : the class which most often calls this routine
    """

    def __init__(self, coords, pot, tol=1e-4, event=None, nsteps=100,
                 nfail_max=200, eigenvec0=None, iprint=-1, orthogZeroEigs=0,
                 nsteps_tangent1=10,
                 nsteps_tangent2=100,
                 lowestEigenvectorQuenchParams=None,
                 tangentSpaceQuenchParams=None,
                 max_uphill_step=0.5,
                 max_uphill_step_initial=0.2,
                 demand_initial_negative_vec=False,
                 negatives_before_check=10,
                 verbosity=1,
                 check_negative=True,
                 invert_gradient=False,
                 hessian_diagonalization=False):
        self.pot = pot
        self.coords = np.copy(coords)
        self.nfev = 0
        self.tol = tol
        self.nsteps = nsteps
        self.event = event
        self.nfail_max = nfail_max
        self.nfail = 0
        self.eigenvec = eigenvec0
        self.orthogZeroEigs = orthogZeroEigs
        self.iprint = iprint
        if lowestEigenvectorQuenchParams is None:
            self.lowestEigenvectorQuenchParams = dict()
        else:
            self.lowestEigenvectorQuenchParams = lowestEigenvectorQuenchParams
        self.max_uphill_step = max_uphill_step
        self.verbosity = verbosity
        self.tangent_space_quencher = mylbfgs  # should make this passable
        if tangentSpaceQuenchParams is None:
            self.tangent_space_quench_params = dict()
        else:
            self.tangent_space_quench_params = dict(tangentSpaceQuenchParams.items())
        self.demand_initial_negative_vec = demand_initial_negative_vec
        self.npositive_max = max(10, self.nsteps // 5)
        self.check_negative = check_negative
        self.negatives_before_check = negatives_before_check
        self.invert_gradient = invert_gradient
        self.hessian_diagonalization = hessian_diagonalization
        if self.verbosity > 0:
            print("will compute the lowest eigenvector by diagonalizing the Hessian")


        self.rmsnorm = 1. / np.sqrt(float(len(coords)))
        self.oldeigenvec = None

        # set tolerance for the tangent space minimization.  
        # Be sure it is at least as tight as self.tol
        self.tol_tangent = self.tol  # * 0.2
        if "tol" in self.tangent_space_quench_params:
            self.tol_tangent = min(self.tol_tangent,
                                   self.tangent_space_quench_params["tol"])
            del self.tangent_space_quench_params["tol"]
        self.tangent_space_quench_params["tol"] = self.tol

        self.nsteps_tangent1 = nsteps_tangent1
        self.nsteps_tangent2 = nsteps_tangent2

        if "maxstep" in self.tangent_space_quench_params:
            self.maxstep_tangent = self.tangent_space_quench_params["maxstep"]
            del self.tangent_space_quench_params["maxstep"]
        else:
            self.maxstep_tangent = 0.1  # this should be determined in a better way

        if "logger" not in self.tangent_space_quench_params:
            self.tangent_space_quench_params["logger"] = logging.getLogger("pele.connect.findTS.tangent_space_quench")


        # set some parameters used in finding lowest eigenvector
        # initial guess for Hermitian
        try:
            self.H0_leig = self.lowestEigenvectorQuenchParams.pop("H0")
        except KeyError:
            self.H0_leig = None

        self.reduce_step = 0
        self.step_factor = .1
        self.npositive = 0

        self._trust_radius = 2.
        self._max_uphill_max = max_uphill_step
        self._max_uphill_min = .01
        if self._max_uphill_min >= self._max_uphill_max:
            self._max_uphill_min = self._max_uphill_max / 5
        self._max_uphill = min(max_uphill_step_initial, self._max_uphill_max)

        self._transverse_walker = None


    @classmethod
    def params(cls, obj=None):
        if obj is None:
            obj = FindTransitionState(np.zeros(2), None)

        params = dict()

        params["tangentSpaceQuenchParams"] = obj.tangent_space_quench_params.copy()
        params["lowestEigenvectorQuenchParams"] = obj.lowestEigenvectorQuenchParams.copy()
        params["tol"] = obj.tol
        params["nsteps"] = obj.nsteps
        params["nfail_max"] = obj.nfail_max
        params["iprint"] = obj.iprint
        params["nsteps_tangent1"] = obj.nsteps_tangent1
        params["nsteps_tangent2"] = obj.nsteps_tangent2
        params["max_uphill_step"] = obj._max_uphill_max
        params["max_uphill_step_initial"] = obj._max_uphill
        params["demand_initial_negative_vec"] = obj.demand_initial_negative_vec
        params["check_negative"] = obj.check_negative
        params["invert_gradient"] = obj.invert_gradient
        params["verbosity"] = obj.verbosity

        # event=None, eigenvec0=None, orthogZeroEigs=0,
        return params

    def _saveState(self, coords):
        """save the state in order to revert a step"""
        self.saved_coords = np.copy(coords)
        self.saved_eigenvec = np.copy(self.eigenvec)
        self.saved_eigenval = self.eigenval
        self.saved_overlap = self.overlap
        self.saved_H0_leig = self.H0_leig
        self.saved_energy = self.energy
        self.saved_gradient = self.gradient.copy()

    # self.saved_oldeigenvec = np.copy(self.oldeigenvec)

    def _resetState(self):
        """restore from a state"""
        coords = np.copy(self.saved_coords)
        self.eigenvec = np.copy(self.saved_eigenvec)
        self.eigenval = self.saved_eigenval
        self.oldeigenvec = np.copy(self.eigenvec)
        self.overlap = self.saved_overlap
        self.H0_leig = self.saved_H0_leig
        self.energy = self.saved_energy
        self.gradient = self.saved_gradient.copy()
        return coords

    def _compute_gradients(self, coords):
        """compute the energy and gradient at the current position and store them for later use"""
        self.nfev += 1
        self.energy, self.gradient = self.pot.getEnergyGradient(coords)

    def get_energy(self):
        """return the already computed energy at the current position"""
        return self.energy

    def get_gradient(self):
        """return the already computed gradient at the current position"""
        return self.gradient

    def run(self):
        """The main loop of the algorithm"""
        coords = np.copy(self.coords)
        res = Result()  # return object
        res.message = []

        self._compute_gradients(coords)
        iend = 0
        for i in range(self.nsteps):
            iend = i
            # get the lowest eigenvalue and eigenvector
            self.overlap = self._getLowestEigenVector(coords, i)
            overlap = self.overlap

            if self.eigenval < 0:
                self.negatives_before_check -= 1

            # determine whether everything looks OK.
            all_ok = self.eigenval < 0 or not self.check_negative
            if not all_ok:
                if i == 0:
                    # we need to accept because we haven't saved the state yet
                    # Also, demand_initial_negative_vec will stop later if needed
                    all_ok = True
            if not all_ok:
                if self.negatives_before_check > 0 and not self.demand_initial_negative_vec:
                    print("  positive before check. setting all ok")
                    all_ok = True

            # if everything is OK, then continue, else revert the step
            if all_ok:
                self._saveState(coords)
                self.reduce_step = 0
            else:
                self.npositive += 1
                if self.npositive > self.npositive_max:
                    logger.warning("positive eigenvalue found too many times. ending %s", self.npositive)
                    res.message.append("positive eigenvalue found too many times %d" % self.npositive)
                    break
                if self.verbosity > 2:
                    logger.info("the eigenvalue turned positive. %s %s", self.eigenval,
                                "Resetting last good values and taking smaller steps")
                coords = self._resetState()
                self.reduce_step += 1

            # step uphill along the direction of the lowest eigenvector
            coords = self._stepUphill(coords)

            # minimize the coordinates in the space perpendicular to the lowest eigenvector
            tangent_ret = self._minimizeTangentSpace(coords, energy=self.get_energy(), gradient=self.get_gradient())
            coords = tangent_ret.coords
            tangentrms = tangent_ret.rms


            # check if we are done and print some stuff
            # self._compute_gradients(coords) # this is unnecessary
            E = self.get_energy()
            grad = self.get_gradient()
            rms = np.linalg.norm(grad) * self.rmsnorm
            gradpar = np.dot(grad, self.eigenvec) / np.linalg.norm(self.eigenvec)

            if self.iprint > 0:
                if (i + 1) % self.iprint == 0:
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
                res.message.append("too many failures in eigenvector search %d" % self.nfail)
                break

            if i == 0 and self.eigenval > 0.:
                if self.verbosity > 1:
                    logger.warning("initial eigenvalue is positive - increase NEB spring constant?")
                if self.demand_initial_negative_vec:
                    logger.warning("            aborting transition state search")
                    res.message.append("initial eigenvalue is positive %f" % self.eigenval)
                    break

        # done.  do one last eigenvector search because coords may have changed
        self._getLowestEigenVector(coords, iend)

        # print some data
        if self.verbosity > 0 or self.iprint > 0:
            logger.info("findTransitionState done: %s %s %s %s %s", iend, E, rms, "eigenvalue", self.eigenval)

        success = True
        # check if results make sense
        if self.eigenval >= 0.:
            if self.verbosity > 2:
                logger.info("warning: transition state is ending with positive eigenvalue %s", self.eigenval)
            success = False
        if rms > self.tol:
            if self.verbosity > 2:
                logger.info("warning: transition state search appears to have failed: rms %s", rms)
            success = False
        if iend >= self.nsteps:
            res.message.append("maximum iterations reached %d" % iend)

        # update nfev with the number of calls from the transverse walker
        if self._transverse_walker is not None:
            twres = self._transverse_walker.get_result()
            self.nfev += twres.nfev

        res.coords = coords
        res.energy = E
        res.eigenval = self.eigenval
        res.eigenvec = self.eigenvec
        res.grad = grad
        res.rms = rms
        res.nsteps = iend
        res.success = success
        res.nfev = self.nfev
        return res

    def _get_lowest_eigenvector_RR(self, coords, gradient=None):
        """get the lowest eigenvector using Reyleigh Ritz minimization"""
        if "nsteps" in self.lowestEigenvectorQuenchParams:
            niter = self.lowestEigenvectorQuenchParams["nsteps"]
        else:
            niter = 100
            if self.verbosity > 3:
                print("Using default of", niter, "steps for finding lowest eigenvalue")
        optimizer = FindLowestEigenVector(coords, self.pot,
                                          eigenvec0=self.eigenvec,
                                          orthogZeroEigs=self.orthogZeroEigs,
                                          gradient=gradient,
                                          **self.lowestEigenvectorQuenchParams)
        # H0=self.H0_leig,
        res = optimizer.run(niter)
        if res.nsteps == 0:
            if self.verbosity > 2:
                print("eigenvector converged, but doing one iteration anyway")
            optimizer.one_iteration()
            res = optimizer.get_result()
        self.H0_leig = res.H0
        return res

    def _get_lowest_eigenvector_diagonalization(self, coords, **kwargs):
        """compute the lowest eigenvector by diagonalizing the Hessian
        
        This scales as N**3, so can be very slow for large systems.
        """
        if self.verbosity > 3:
            print("computing the lowest eigenvector by diagonalizing the Hessian")
        hess = self.pot.getHessian(coords)
        eigenval, evec = get_smallest_eig(hess)
        res = Result()
        res.eigenval = eigenval
        res.eigenvec = evec
        res.nfev = 1
        res.success = True
        res.rms = 0.
        return res


    def _getLowestEigenVector(self, coords, i, gradient=None):
        """compute the lowest eigenvector at position coords
        
        Parameters
        ----------
        coords : the current position
        i : the iteration number
        gradient : the gradient at coords
        """
        if self.hessian_diagonalization:
            res = self._get_lowest_eigenvector_diagonalization(coords)
        else:
            res = self._get_lowest_eigenvector_RR(coords, gradient=gradient)

        self.leig_result = res
        self.nfev += res.nfev

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
        """minimize the energy in the space perpendicular to eigenvec.
        
        Parameters
        ----------
        coords : the current position
        energy, gradient : the energy and gradient at the current position
        """
        assert gradient is not None
        if self._transverse_walker is None:
            if self.invert_gradient:
                # note: if we pass transverse energy and gradient here we can save 1 potential call
                self._transverse_walker = _DimerTranslator(coords, self.pot, self.eigenvec,
                                                           **self.tangent_space_quench_params)
            else:
                self._transverse_walker = _TransverseWalker(coords, self.pot, self.eigenvec, energy, gradient,
                                                            **self.tangent_space_quench_params)
        else:
            self._transverse_walker.update_eigenvec(self.eigenvec, self.eigenval)
            self._transverse_walker.update_coords(coords, energy, gradient)

        # determine the number of steps
        # i.e. if the eigenvector is deemed to have converged or is changing slowly
        eigenvec_converged = self.overlap > .999
        if eigenvec_converged:
            nstepsperp = self.nsteps_tangent2
        else:
            nstepsperp = self.nsteps_tangent1

        # reduce the maximum step size if necessary
        maxstep = self.maxstep_tangent
        if self.reduce_step > 0:
            maxstep *= self.step_factor ** self.reduce_step
        self._transverse_walker.update_maxstep(maxstep)

        coords_old = coords.copy()
        ret = self._transverse_walker.run(nstepsperp)

        coords = ret.coords
        self.tangent_move_step = np.linalg.norm(coords - coords_old)
        self.tangent_result = ret
        if self.tangent_move_step > 1e-16:
            try:
                self.energy, self.gradient = self._transverse_walker.get_true_energy_gradient(coords)
            except AttributeError:
                raise Exception("was tspot was never called? use the same gradient")
        return ret

    def _update_max_uphill_step(self, Fold, stepsize):
        """use a trust radius to update the maximum uphill step size"""
        Fnew = np.dot(self.eigenvec, self.get_gradient())
        # EPER=MIN(DABS(1.0D0-(FOBNEW-FOB)/(PSTEP*EVALMIN)),DABS(1.0D0-(-FOBNEW-FOB)/(PSTEP*EVALMIN)))
        a1 = 1. - (Fnew - Fold) / (stepsize * self.eigenval)
        a2 = 1. - (-Fnew - Fold) / (stepsize * self.eigenval)
        eper = min(np.abs(a1), np.abs(a2))
        if eper > self._trust_radius:
            # reduce the maximum step size
            self._max_uphill = max(self._max_uphill / 1.1, self._max_uphill_min)
            if self.verbosity > 2:
                print("decreasing max uphill step to", self._max_uphill, "Fold", Fold, "Fnew", Fnew, "eper", eper, "eval", self.eigenval)
        else:
            # increase the maximum step size
            self._max_uphill = min(self._max_uphill * 1.1, self._max_uphill_max)
            if self.verbosity > 2:
                print("increasing max uphill step to", self._max_uphill, "Fold", Fold, "Fnew", Fnew, "eper", eper, "eval", self.eigenval)


    def _stepUphill(self, coords):
        """step uphill in the direction of self.eigenvec.
        
        this is an eigenvector following step uphill.  The equation for the step size
        is described in 
        
        DJ Wales, J chem phys, 1994, 101, 3750--3762
        Rearrangements of 55-atom lennard-jones and (c-60)(55) clusters
        
        self.eigenval is used to determine the best stepsize
        """
        # the energy and gradient are already known
        grad = self.get_gradient()
        F = np.dot(grad, self.eigenvec)
        h = 2. * F / np.abs(self.eigenval) / (1. + np.sqrt(1. + 4. * (F / self.eigenval) ** 2))

        if self.eigenval > 0 and self.verbosity >= 2:
            logger.warn("eigenvalue is positive, but stepping uphill along the lowest curvature mode anyway")

        # get the maxstep and scale it if necessary
        maxstep = self._max_uphill
        if self.reduce_step > 0:
            maxstep *= self.step_factor ** self.reduce_step

        if np.abs(h) > maxstep:
            if self.verbosity >= 5:
                logger.debug("reducing uphill step from %s %s %s", h, "to", maxstep)
            h *= maxstep / np.abs(h)
        self.uphill_step_size = h
        coords += h * self.eigenvec

        # recompute the energy and gradient
        Eold = self.energy
        self._compute_gradients(coords)

        if self.energy < Eold and self.verbosity > 0:
            logger.warn("energy decreased after uphill step %s -> %s", Eold, self.energy)


        # update the maximum step using a trust ratio
        if self.eigenval < 0:
            self._update_max_uphill_step(F, h)

        if self.verbosity > 2:
            logger.info("stepping uphill with stepsize %s", h)

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

