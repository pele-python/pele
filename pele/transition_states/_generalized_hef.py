"""
routines for a generalized hybrid eigenvector following
"""
from __future__ import print_function
import numpy as np

from pele.transition_states._transverse_walker import _TransversePotential
from pele.optimize import MYLBFGS, Result


class _HybridEigenvectorWalker(object):
    """a class to perform the translational steps in hybrid eigenvector following
    
    Notes
    -----
    step uphill along the direction of least curvature (smallest eigenvector) using
    Newton's method
    
    minimize in the space perpendicular to the smallest eigenvector using
    an optimization algorithm
    """

    def __init__(self, coords, potential, eigenvec, eigenval=None,
                 max_uphill_step=0.1,
                 nsteps_tangent=10,
                 verbosity=5,
                 tol=1e-4,
                 **transverse_kwargs
    ):
        print("in HEF")
        self.eigenval = eigenval
        self.nsteps_tangent = nsteps_tangent
        self.max_uphill_step = max_uphill_step
        self.coords = coords
        self.verbosity = verbosity
        self.transverse_kwargs = transverse_kwargs
        self.tol = tol
        if self.transverse_kwargs is None: self.transverse_kwargs = dict()
        if "tol" not in self.transverse_kwargs:
            self.transverse_kwargs["tol"] = self.tol / 2.

        self.transverse_potential = _TransversePotential(potential, eigenvec)
        self.update_eigenvec(eigenvec, eigenval)
        # self.transverse_walker = MYLBFGS(self.coords, self.potential, **transverse_kwargs)
        self._transverse_walker_state = None

        # self.energy, self.grad = self.potential.getEnergyGradient(self.coords)
        self._update_coords(coords)

        self.iter_number = 0
        self.debug = True

    def stop_criterion_satisfied(self):
        rms = np.linalg.norm(self.gradient) / np.sqrt(self.coords.size)
        return rms < self.tol


    def _update_coords(self, coords, transverse_energy=None,
                       transverse_gradient=None):
        self.coords = coords
        if transverse_energy is None or transverse_gradient is None:
            self._transverse_energy, self._transverse_gradient = self.transverse_potential.getEnergyGradient(coords)
        else:
            self._transverse_energy = transverse_energy
            self._transverse_gradient = transverse_gradient
        self.energy = self.transverse_potential._true_energy
        self.gradient = self.transverse_potential._true_gradient.copy()

    def update_eigenvec(self, eigenvec, eigenval):
        self.transverse_potential.update_vector(eigenvec)
        self.eigenval = eigenval

    def get_gradient(self):
        return self.gradient

    def get_energy(self):
        return self.energy

    def get_eigenvector(self):
        return self.transverse_potential.vector

    def _step_uphill(self, coords, gradient, eigenvec, eigenval):
        """step uphill along the direction of self.eigenvec
        
        """
        F = np.dot(gradient, eigenvec)
        stepsize = 2. * F / np.abs(self.eigenval) / (1. + np.sqrt(1. + 4. * F ** 2 / eigenval ** 2))

        maxstep = self.max_uphill_step
        # if self.reduce_step > 0:
        # maxstep *= (self.step_factor)**self.reduce_step


        if np.abs(stepsize) > maxstep:
            if self.verbosity >= 5:
                # logger.debug("reducing step from %s %s %s", stepsize, "to", maxstep)
                print("reducing uphill step from %s %s %s" % (np.abs(stepsize), "to", maxstep))
            stepsize *= maxstep / np.abs(stepsize)

        print("uphill step", stepsize)
        if self.debug:
            ne, ng = self.transverse_potential.pot.getEnergyGradient(coords + stepsize * eigenvec)
            print("uphill step", stepsize, "dE", ne - self.energy, "gradpar old new", \
                np.dot(self.gradient, self.get_eigenvector()), np.dot(ng, self.get_eigenvector()))

        return coords + stepsize * eigenvec


    def _minimize_transverse(self, nsteps, transverse_energy=None, transverse_gradient=None):
        # must update
        minimizer = MYLBFGS(self.coords, self.transverse_potential,
                            nsteps=self.nsteps_tangent,
                            energy=transverse_energy, gradient=transverse_gradient,
                            **self.transverse_kwargs)

        if self._transverse_walker_state is not None:
            minimizer.set_state(self._transverse_walker_state)
        ret = minimizer.run()
        self._transverse_walker_state = minimizer.get_state()
        self._stop_criterion_satisfied = minimizer.stop_criterion_satisfied()

        return ret


    def _print_status(self):
        rmsnorm = 1. / np.sqrt(self.coords.size)
        rmstv = np.linalg.norm(self._transverse_gradient) * rmsnorm
        rmstrue = np.linalg.norm(self.gradient) * rmsnorm
        print("rms true", rmstrue, "transverse", rmstv, "gradpar", np.dot(self.gradient, self.get_eigenvector()), \
            "eigenval", self.eigenval)


    def one_iteration(self):
        coords = self._step_uphill(self.coords, self.gradient,
                                   self.get_eigenvector(), self.eigenval)
        self._update_coords(coords)
        ret = self._minimize_transverse(self.nsteps_tangent, transverse_energy=self._transverse_energy,
                                        transverse_gradient=self._transverse_gradient)

        # update some values for the next iteration
        self._update_coords(ret.coords.copy(), transverse_energy=ret.energy,
                            transverse_gradient=ret.grad)
        # self.coords = ret.coords.copy()
        # self._transverse_energy = ret.energy
        # self._transverse_gradient = ret.grad
        # self.energy = self.transverse_potential.true_energy
        # self.gradient = self.transverse_potential.true_gradient

        self._print_status()

    def run(self, niter):
        self.nsteps_tangent = niter
        self.one_iteration()

    def get_result(self):
        res = Result()
        res.energy = self.energy
        res.gradient = self.gradient
        res.coords = self.coords
        res.nsteps = self.iter_number
        res.rms = np.linalg.norm(res.gradient) / np.sqrt(len(res.gradient))
        res.nfev = self.transverse_potential.nfev

        res.eigenval = self.eigenval
        res.eigenvec = self.get_eigenvector()

        res.success = self.stop_criterion_satisfied()
        return res

    def get_true_energy_gradient(self, coords):
        # js850> I added this function without carefully checking it
        return self.transverse_potential.get_true_energy_gradient(coords)


#
# testing only
#

def test():  # pragma: no cover
    from pele.systems import LJCluster
    from pele.transition_states import GeneralizedDimer
    from pele.transition_states import FindTransitionState
    from pele.utils.xyz import read_xyz

    system = LJCluster(13)
    x = system.get_random_configuration()
    x = read_xyz(open("tests/lj18_ts.xyz", "r")).coords.flatten()

    dimer = GeneralizedDimer(x.copy(), system.get_potential(),
                             rotational_steps=40,
                             dimer=False,
                             leig_kwargs=dict(iprint=10, maxstep=2.,
                                              tol=.1,
                             ),
                             translator_kwargs=dict(iprint=1,
                                                    verbosity=10,
                             ),
    )
    ret = dimer.run()
    print(ret)

    print("\n\nnow the same with the old version")
    searcher = FindTransitionState(x.copy(), system.get_potential(), iprint=1, verbosity=10)
    # print ret
    ret = searcher.run()
    print(ret)


if __name__ == "__main__":
    test()

