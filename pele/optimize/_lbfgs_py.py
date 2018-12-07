from __future__ import print_function
from __future__ import absolute_import
import numpy as np
import logging
from collections import namedtuple

from .optimization_exceptions import LineSearchError
from pele.optimize import Result

__all__ = ["LBFGS"]

_logger = logging.getLogger("pele.optimize")


class LBFGS(object):
    """
    minimize a function using the LBFGS routine
    
    Parameters
    ----------
    X : array
        the starting configuration for the minimization
    pot :
        the potential object
    nsteps : int
        the maximum number of iterations
    tol : float
        the minimization will stop when the rms grad is less than tol
    iprint : int
        how often to print status information
    maxstep : float
        the maximum step size
    maxErise : float
        the maximum the energy is alowed to rise during a step.
        The step size will be reduced until this condition is satisfied.
    M : int
        the number of previous iterations to use in determining the optimal step
    rel_energy : bool
        if True, then maxErise the the *relative* maximum the energy is allowed
        to rise during a step
    H0 : float
        the initial guess for the inverse diagonal Hessian.  This particular
        implementation of LBFGS takes all the inverse diagonal components to be the same. 
    events : list of callables
        these are called after each iteration.  events can also be added using
        attachEvent()
    alternate_stop_criterion : callable
        this criterion will be used rather than rms gradiant to determine when
        to stop the iteration
    debug : 
        print debugging information
    logger : logger object
        messages will be passed to this logger rather than the default
    energy, gradient : float, float array
        The initial energy and gradient.  If these are both not None then the
        energy and gradient of the initial point will not be calculated, saving
        one potential call.
    armijo : bool
        Use the Armijo criterion instead of maxErise as the criterion for 
        accepting a step size.  The Armijo criterion is the first wolfe criterion
        and is a condition that the energy decrease sufficiently
    armijo_c : float
        This adjusts how strong the armijo rule is.  0 < armijo_c < 1.  Default 
        1e-4
    fortran : bool
        use the fortran version of the LBFGS.  Only the step which computes
        the step size and direction from the memory is in fortran.
         
    Notes
    -----
    This each iteration of this minimization routine is composed of the 
    following parts
    
    1. determine a step size and direction using the LBFGS algorithm
    
    2. ensure the step size is appropriate (see maxErise and maxstep).
       Reduce the step size until conditions are satisfied.
    
    3. take step

    http://dx.doi.org/10.1007/BF01589116
    
    See Also
    --------
    MYLBFGS : this implemented in a compiled language
    lbfgs_py : a function wrapper
    
    """

    def __init__(self, X, pot, maxstep=0.1, maxErise=1e-4, M=4,
                 rel_energy=False, H0=0.1, events=None,
                 alternate_stop_criterion=None, debug=False,
                 iprint=-1, nsteps=10000, tol=1e-5, logger=None,
                 energy=None, gradient=None, armijo=False,
                 armijo_c=1e-4,
                 fortran=False):
        X = X.copy()
        self.X = X
        self.N = len(X)
        self.M = M
        self.pot = pot
        self._use_wolfe = False  # this didn't work very well.  should probably remove
        self._armijo = bool(armijo)
        self._wolfe1 = armijo_c
        self._wolfe2 = 0.99
        self._cython = False  # we could make this passable
        self._fortran = bool(fortran)
        self.funcalls = 0
        if energy is not None and gradient is not None:
            self.energy = energy
            self.G = gradient
        else:
            self.energy, self.G = self.pot.getEnergyGradient(self.X)
            self.funcalls += 1
        self.rms = np.linalg.norm(self.G) / np.sqrt(self.N)

        self.maxstep = maxstep
        self.maxErise = maxErise
        self.rel_energy = rel_energy  # use relative energy comparison for maxErise
        self.events = events  # a list of events to run during the optimization
        if self.events is None: self.events = []
        self.iprint = iprint
        self.nsteps = nsteps
        self.tol = tol
        if logger is None:
            self.logger = _logger
        else:
            self.logger = logger

        self.alternate_stop_criterion = alternate_stop_criterion
        self.debug = debug  # print debug messages

        self.s = np.zeros([self.M, self.N])  # position updates
        self.y = np.zeros([self.M, self.N])  # gradient updates

        if H0 is None:
            self.H0 = 0.1
        else:
            self.H0 = H0
        if self.H0 < 1e-10:
            self.logger.warning("initial guess for inverse Hessian diagonal is negative or too small %s %s",
                                self.H0, "resetting it to 0.1")
            self.H0 = 0.1
        self.rho = np.zeros(M)
        self.k = 0

        self.s[0, :] = self.X
        self.y[0, :] = self.G
        self.rho[0] = 0.  # 1. / np.dot(X,G)

        self.dXold = np.zeros(self.X.size)
        self.dGold = np.zeros(self.X.size)
        self._have_dXold = False

        self.nfailed = 0

        self.iter_number = 0
        self.result = Result()
        self.result.message = []

    def get_state(self):
        """return the state of the LBFGS memory"""
        State = namedtuple("State", "s y rho k H0 dXold dGold have_dXold")
        state = State(s=self.s.copy(), y=self.y.copy(),
                      rho=self.rho.copy(), k=self.k, H0=self.H0,
                      dXold=self.dXold.copy(), dGold=self.dGold.copy(),
                      have_dXold=self._have_dXold)
        return state

    def set_state(self, state):
        """set the LBFGS memory from the passed state"""
        self.s = state.s
        self.y = state.y
        self.rho = state.rho
        self.k = state.k
        self.H0 = state.H0
        self.dXold = state.dXold
        self.dGold = state.dGold
        self._have_dXold = state.have_dXold
        assert self.s.shape == (self.M, self.N)
        assert self.y.shape == (self.M, self.N)
        assert self.rho.shape == (self.M,)
        assert self.dXold.shape == (self.N,)
        assert self.dGold.shape == (self.N,)

    def update_coords(self, X, E, G):
        """change the location of the minimizer manually
        
        If the position change is too great the LBFGS memory will not accurately
        represent the local curvature. 
        """
        self.X = X.copy()
        self.energy = float(E)
        self.G = G.copy()
        self.rms = np.linalg.norm(self.G) / np.sqrt(self.G.size)

    def _add_step_to_memory(self, dX, dG):
        """
        add a step to the LBFGS memory
        
        Parameters
        ----------
        dX : ndarray
            the step: X - Xold
        dG : ndarray 
            the change in gradient along the step: G - Gold
        """
        klocal = (self.k + self.M) % self.M  # =k  cyclical
        self.s[klocal, :] = dX
        self.y[klocal, :] = dG

        # the local curvature along direction s is np.dot(y, s) / norm(s)
        YS = np.dot(dX, dG)
        if YS == 0.:
            self.logger.warning("resetting YS to 1 in lbfgs %s", YS)
            YS = 1.
        self.rho[klocal] = 1. / YS

        # update the approximation for the diagonal inverse hessian
        # scale H0 according to
        # H_k = YS/YY * H_0 
        # this is described in Liu and Nocedal 1989 
        # http://dx.doi.org/10.1007/BF01589116
        # note: for this step we assume H0 is always the identity
        # js850: This ability to update H0 is what sets LBFGS apart from BFGS
        # and makes it such a superior algorithm in my opinion.  This is why
        # LBFGS gets away with not using a more complicated linesearch algorithm
        # and why BFGS (which can't have this step) gives nonsensical results without
        # a linesearch.
        YY = np.dot(dG, dG)
        if YY == 0.:
            self.logger.warning("warning: resetting YY to 1 in lbfgs %s", YY)
            YY = 1.
        self.H0 = YS / YY

        # increment k
        self.k += 1

    def _get_LBFGS_step_cython(self, G):
        from . import _cython_lbfgs

        return _cython_lbfgs._compute_LBFGS_step(G, self.s, self.y, self.rho,
                                                 self.k, self.H0)

    def _get_LBFGS_step_fortran(self, G):
        from . import mylbfgs_updatestep

        ret = mylbfgs_updatestep.lbfgs_get_step_wrapper(G, self.s.reshape(-1), self.y.reshape(-1), self.rho,
                                                        self.k, self.H0)
        return ret

    def _get_LBFGS_step(self, G):
        """use the LBFGS algorithm to compute a suggested step from the memory
        """
        if self._cython:
            return self._get_LBFGS_step_cython(G)
        elif self._fortran:
            return self._get_LBFGS_step_fortran(G)
        s = self.s
        y = self.y
        rho = self.rho
        k = self.k

        q = G.copy()
        a = np.zeros(self.M)
        myrange = [i % self.M for i in range(max([0, k - self.M]), k, 1)]
        assert len(myrange) == min(self.M, k)
        for i in reversed(myrange):
            a[i] = rho[i] * np.dot(s[i, :], q)
            q -= a[i] * y[i, :]

        # z[:] = self.H0[ki] * q[:]
        z = q  # q is not used anymore after this, so we can use it as workspace
        z *= self.H0
        for i in myrange:
            beta = rho[i] * np.dot(y[i, :], z)
            z += s[i, :] * (a[i] - beta)

        stp = -z

        if k == 0:
            # make first guess for the step length cautious
            gnorm = np.linalg.norm(G)
            stp *= min(gnorm, 1. / gnorm)

        return stp

    def getStep(self, X, G):
        """update the LBFGS memory and compute a step direction and size
        
        http://en.wikipedia.org/wiki/Limited-memory_BFGS
        
        Liu and Nocedal 1989
        http://dx.doi.org/10.1007/BF01589116
        """
        # we have a new X and G, save in s and y
        if self._have_dXold:
            self._add_step_to_memory(self.dXold, self.dGold)

        stp = self._get_LBFGS_step(G)

        return stp

    def adjustStepSize(self, X, E, G, stp):
        """
        We now have a proposed step.  This function will make sure it is 
        a good step and then take it.  This is known as a Backtracking linesearch
        
        http://en.wikipedia.org/wiki/Backtracking_line_search
        
        1) if the step is not anti-aligned with the gradient (i.e. downhill), 
           then reverse the step
        
        2) if the step is larger than maxstep, then rescale the step
        
        3) calculate the energy and gradient of the new position
        
        4) if the step increases the energy by more than maxErise, 
            then reduce the step size and go to 3)
        
        5) if the step is reduced more than 10 times and the energy is still 
           not acceptable, then increment nfail, reset the lbfgs optimizer and 
           continue
            
        6) if nfail is greater than 5 abort the quench
                
        """
        f = 1.
        X0 = X.copy()
        G0 = G.copy()
        E0 = E

        if np.dot(G, stp) > 0:
            if self.debug:
                overlap = np.dot(G, stp) / np.linalg.norm(G) / np.linalg.norm(stp)
                self.logger.warn("LBFGS returned uphill step, reversing step direction: overlap {}".format(overlap))
            stp = -stp

        stepsize = np.linalg.norm(stp)

        if f * stepsize > self.maxstep:
            f = self.maxstep / stepsize

        nincrease = 0
        while True:
            X = X0 + f * stp
            E, G = self.pot.getEnergyGradient(X)
            self.funcalls += 1

            # if the increase is greater than maxErise reduce the step size
            if self._accept_step(E, E0, G, G0, f * stp):
                break
            else:
                if self.debug:
                    self.logger.warn("energy increased, trying a smaller step %s %s %s %s", E, E0, f * stepsize,
                                     nincrease)
                f /= 10.
                nincrease += 1
                if nincrease > 10:
                    break

        if nincrease > 10:
            self.nfailed += 1
            if self.nfailed > 10:
                raise LineSearchError("lbfgs: too many failures in adjustStepSize, exiting")

            # abort the linesearch, reset the memory and reset the coordinates            
            self.logger.warning("lbfgs: having trouble finding a good step size. %s %s, resetting the minimizer",
                                f * stepsize, stepsize)
            self.reset()
            E = E0
            G = G0
            X = X0
            f = 0.

        self.stepsize = f * stepsize
        return X, E, G

    def _accept_step(self, Enew, Eold, Gnew, Gold, step, strong=False):
        """determine whether the step is acceptable"""
        if self._use_wolfe:
            return self._wolfe_conditions(Enew, Eold, Gnew, Gold, step, strong)
        elif self._armijo:
            return self._armijo_condition(Enew, Eold, Gold, step)
        else:
            # get the increase in energy            
            if self.rel_energy:
                if Eold == 0:
                    Eold = 1e-100
                dE = (Enew - Eold) / abs(Eold)
            else:
                dE = Enew - Eold

            # if the increase is greater than maxErise reduce the step size
            return dE <= self.maxErise

    def _armijo_condition(self, Enew, Eold, Gold, step, return_overlap=False):
        """test if the armijo condition is satisfied
        
        The energy cannot rise more than an amount dependent on the 
        dot product of the gradient and the step
        """
        overlap_old = np.dot(Gold, step)
        armijo = Enew <= Eold + overlap_old * self._wolfe1
        if not armijo and self.debug:
            stepsize = np.linalg.norm(step)
            print(self.iter_number, "rejecting step due to energy", Enew, Enew - Eold, overlap_old * self._wolfe1, "stepsize", stepsize)
        if return_overlap:
            return armijo, overlap_old
        return armijo

    def _wolfe_conditions(self, Enew, Eold, Gnew, Gold, step, strong=False):
        """return True if the Wolfe conditions are satisfied, False otherwise
        
        wolfe1 : the energy cannot rise more than an amount dependent on the 
                 dot product of the gradient and the step
        
        wolfe2 : the overlap of the gradient with the step direction cannot 
        decrease by more than a given factor 
        """
        armijo, overlap_old = self._armijo_condition(Enew, Eold, Gold, step, return_overlap=True)
        if not armijo:
            return False

        overlap_new = np.dot(Gnew, step)
        if strong:
            wolfe2 = np.abs(overlap_new) <= np.abs(overlap_old) * self._wolfe2
        else:
            wolfe2 = overlap_new >= overlap_old * self._wolfe2
        if not wolfe2 and self.debug:
            stepsize = np.linalg.norm(step)
            print("wolfe:", self.iter_number, "rejecting step due to gradient", overlap_new, overlap_old, self._wolfe2, "stepsize", stepsize)
        return armijo and wolfe2


    def reset(self):
        """reset the LBFGS memory and H0"""
        self.H0 = 0.1
        self.k = 0
        self._have_dXold = False

    def attachEvent(self, event):
        self.events.append(event)

    def one_iteration(self):
        """do one iteration of the LBFGS loop
        """
        stp = self.getStep(self.X, self.G)

        Xnew, self.energy, Gnew = self.adjustStepSize(self.X, self.energy, self.G, stp)
        self.dXold = Xnew - self.X
        self.dGold = Gnew - self.G
        self._have_dXold = True
        self.X = Xnew
        self.G = Gnew

        self.rms = np.linalg.norm(self.G) / np.sqrt(self.N)

        if self.iprint > 0 and self.iter_number % self.iprint == 0:
            self.logger.info("lbfgs: %s %s %s %s %s %s %s %s %s", self.iter_number, "E", self.energy,
                             "rms", self.rms, "funcalls", self.funcalls, "stepsize", self.stepsize)
        for event in self.events:
            event(coords=self.X, energy=self.energy, rms=self.rms)

        self.iter_number += 1
        return True

    def stop_criterion_satisfied(self):
        """test the stop criterion"""
        if self.alternate_stop_criterion is None:
            return self.rms < self.tol
        else:
            return self.alternate_stop_criterion(energy=self.energy, gradient=self.G,
                                                 tol=self.tol, coords=self.X)

    def run(self):
        """run the LBFGS minimizer
        
        stop when the stop criterion is satisfied or  when the maximum number 
        of steps is reached
        
        Returns
        -------
        return a results object
        """
        while self.iter_number < self.nsteps and not self.stop_criterion_satisfied():
            try:
                self.one_iteration()
            except LineSearchError:
                self.logger.error("problem with adjustStepSize, ending quench")
                self.rms = np.linalg.norm(self.G) / np.sqrt(self.N)
                self.logger.error("    on failure: quench step %s %s %s %s", self.iter_number, self.energy, self.rms,
                                  self.funcalls)
                self.result.message.append("problem with adjustStepSize")
                break

        return self.get_result()

    def get_result(self):
        """return a results object"""
        res = self.result
        res.nsteps = self.iter_number
        res.nfev = self.funcalls
        res.coords = self.X
        res.energy = self.energy
        res.rms = self.rms
        res.grad = self.G
        res.H0 = self.H0
        res.success = self.stop_criterion_satisfied()
        return res

#
# only testing stuff below here
#
#
#class PrintEvent:
#    def __init__(self, fname):
#        self.fout = open(fname, "w")
#        self.coordslist = []
#
#    def __call__(self, coords, **kwargs):
#        from pele.utils.xyz import write_xyz
#        write_xyz(self.fout, coords)
#        self.coordslist.append( coords.copy() )
#        
#    
#
#def test(pot, natoms = 100, iprint=-1):    
#    #X = bfgs.getInitialCoords(natoms, pot)
#    #X += np.random.uniform(-1,1,[3*natoms]) * 0.3
#    X = np.random.uniform(-1,1,[natoms*3])*(1.*natoms)**(1./3)*.5
#    
#    runtest(X, pot, natoms, iprint)
#
#def runtest(X, pot, natoms = 100, iprint=-1):
#    tol = 1e-5
#
#    Xinit = np.copy(X)
#    e, g = pot.getEnergyGradient(X)
#    print "energy", e
#    
#    lbfgs = LBFGS(X, pot, maxstep = 0.1, tol=tol, iprint=iprint, nsteps=10000)
#    printevent = PrintEvent( "debugout.xyz")
#    lbfgs.attachEvent(printevent)
#    
#    ret = lbfgs.run()
#    print "done", ret
#    
#    print "now do the same with scipy lbfgs"
#    from pele.optimize import lbfgs_scipy as quench
#    ret = quench(Xinit, pot, tol = tol)
#    print ret 
#    
#    if False:
#        print "now do the same with scipy bfgs"
#        from pele.optimize import bfgs as oldbfgs
#        ret = oldbfgs(Xinit, pot, tol = tol)
#        print ret    
#    
#    if False:
#        print "now do the same with gradient + linesearch"
#        import _bfgs
#        gpl = _bfgs.GradientPlusLinesearch(Xinit, pot, maxstep = 0.1)  
#        ret = gpl.run(1000, tol = 1e-6)
#        print ret 
#            
#    if False:
#        print "calling from wrapper function"
#        from pele.optimize import lbfgs_py
#        ret = lbfgs_py(Xinit, pot, tol = tol)
#        print ret
#
#
#    try:
#        import pele.utils.pymolwrapper as pym
#        pym.start()
#        for n, coords in enumerate(printevent.coordslist):
#            coords = coords.reshape([-1, 3])
#            pym.draw_spheres(coords, "A", n)
#    except ImportError:
#        print "error loading pymol"
#
#        
#if __name__ == "__main__":
#    from pele.potentials.lj import LJ
#    from pele.potentials.ATLJ import ATLJ
#    pot = ATLJ()
#
#    #test(pot, natoms=3, iprint=1)
#    
##    coords = np.loadtxt("coords")
#    natoms = 10
#    coords = np.random.uniform(-1,1,natoms*3)
#    print coords.size
#    coords = np.reshape(coords, coords.size)
#    print coords
#    runtest(coords, pot, natoms=3, iprint=1)
#    
#
#









