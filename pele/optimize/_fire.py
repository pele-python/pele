"""
Created on 30 Apr 2012

@author: ruehle
"""
from __future__ import print_function
import numpy as np
import math
import logging

from pele.optimize import Result

__all__ = ["Fire"]

_logger = logging.getLogger("pele.optimize")


class Fire(object):
    """
    The FIRE optimization algorithm

    Parameters
    ----------
    coords : array
        the starting configuration
    potential :
        the potential

    Parameters
    ----------
    coords : array
        the starting configuration for the optimization.
    potential :
        the potential object
    dt : float
        adaptive time step for integration of the equations of motion
    maxstep : float
        the maximum step size permitted
    dtmax : float
        the maximum time step permitted
    Nmin : int
        minimum number of steps taken after an uphill step before
        adaptive parameters are allowed to change
    finc : float
        factor by which the time step is increased if moving downhill and
        at least Nmin steps from the last uphill step
    fdec : float
        factor by which the time step is decreased if an uphill step has been taken
    astart : float
        initial value of a and value taken after an uphill step.
    fa : float
        factor by which parameter a decreases if moving downhill and at
        least Nmin steps from the last uphill step
    a : float
        adaptive parameter that controls the velocity used to evolve the system.
    iprint : int
        how often to print status information
    alternative_stop_criterion : callable
        this criterion will be used rather than rms gradient to determine
        when to stop the iteration.  must return a boolean and accept the following keywords:

            1. `energy`
            2. `gradient`
            3. `tol`
    events : list of callables
        these are called after each iteration. Events can also be added using Fire.attachEvent().
        Each event must accept keywords:

            1. `coords`
            2. `energy`
            3. `rms`

    Notes
    -----
    The Fast Inertial Relaxation Engine is an optimization algorithm based
    on molecular dynamics with modifications to the velocity and adaptive
    time steps. The method is based on a blind skier searching for the
    bottom of a valley and is described and tested here:

    Erik Bitzek, Pekka Koskinen, Franz Gaehler, Michael Moseler, and Peter Gumbsch.
    Phys. Rev. Lett. 97, 170201 (2006)
    http://link.aps.org/doi/10.1103/PhysRevLett.97.170201

    This implementation of the algorithm differs significantly from the original
    algorithm in the order in which the steps are taken.

    See Also
    --------
    LBFGS
    MYLBFGS
    """

    def __init__(self, coords, potential, restart=None, logfile='-', trajectory=None,
                 dt=0.1, maxstep=0.5, dtmax=1., Nmin=5, finc=1.1, fdec=0.5,
                 astart=0.1, fa=0.99, a=0.1, iprint=-1,
                 alternate_stop_criterion=None, events=None, logger=None):

        self.dt = dt
        self.Nsteps = 0
        self.maxstep = maxstep
        self.dtmax = dtmax
        self.Nmin = Nmin
        self.finc = finc
        self.fdec = fdec
        self.astart = astart
        self.fa = fa
        self.a = a
        self.coords = coords
        self.potential = potential
        self.v = None
        self.nsteps = 0
        self.iprint = iprint
        self.alternate_stop_criterion = alternate_stop_criterion
        if logger is None:
            self.logger = _logger
        else:
            self.logger = logger
        if events is None:
            self.events = []
        else:
            self.events = events

        self.nfev = 0

    def initialize(self):
        self.v = None

    def step(self, f):
        coords = self.coords
        if self.v is None:
            self.v = np.zeros((len(coords)))
        else:
            vf = np.vdot(f, self.v)
            if vf > 0.0:
                self.v = (1.0 - self.a) * self.v + self.a * f / np.sqrt(
                    np.vdot(f, f)) * np.sqrt(np.vdot(self.v, self.v))
                if self.Nsteps > self.Nmin:
                    self.dt = min(self.dt * self.finc, self.dtmax)
                    self.a *= self.fa
                self.Nsteps += 1
            else:
                self.v[:] *= 0.0
                self.a = self.astart
                self.dt *= self.fdec
                self.Nsteps = 0

        self.v += self.dt * f
        dr = self.dt * self.v
        if False:  # how do we determine maxstep?
            normdr = np.sqrt(np.vdot(dr, dr))
        else:
            normdr = max(np.abs(dr))
        if normdr > self.maxstep:
            dr = self.maxstep * dr / normdr
        self.coords = coords + dr

    def run(self, fmax=1e-3, steps=100000):
        """Run structure optimization algorithm.

        This method will return when the forces on all individual
        atoms are less than *fmax* or when the number of steps exceeds
        *steps*.
        """
        self.fmax = fmax
        step = 0
        res = Result()
        res.success = False
        while step < steps:
            E, f = self.potential.getEnergyGradient(self.coords)
            self.nfev += 1
            if self.alternate_stop_criterion is None:
                i_am_done = self.converged(f)
            else:
                i_am_done = self.alternate_stop_criterion(energy=E, gradient=f,
                                                          tol=self.fmax)
            if i_am_done:
                res.success = True
                break
            self.step(-f)
            self.nsteps += 1
            rms = np.linalg.norm(f) / np.sqrt(len(f))
            if self.iprint > 0:
                if step % self.iprint == 0:
                    self.logger.info("fire: %s E %s rms %s", step, E, rms)
            for event in self.events:
                event(coords=self.coords, energy=E, rms=rms)

            step += 1

        res.nsteps = step
        res.nfev = self.nfev
        res.coords = self.coords
        res.energy = E
        res.grad = -f
        res.rms = np.linalg.norm(res.grad) / np.sqrt(len(res.grad))
        self.result = res
        return res


    def converged(self, forces=None):
        """Did the optimization converge?"""
        return np.linalg.norm(forces) / math.sqrt(len(forces)) < self.fmax


if __name__ == "__main__":
    import pele.potentials.lj as lj

    pot = lj.LJ()
    coords = 10. * np.random.random(300)
    opt = Fire(coords, pot.getEnergyGradient, dtmax=0.1, dt=0.01, maxstep=0.01, iprint=200)
    opt.run(fmax=1e-1, steps=10000)
    print(pot.getEnergy(opt.coords))
    print(opt.nsteps)

