from __future__ import print_function
__all__ = ["AdaptiveStepsize"]

from .generic import TakestepInterface


class AdaptiveStepsize(TakestepInterface):
    """Adaptive stepsize adjustment

    AdaptiveStepsize is a wrapper for a takestep object and automatically adjusts
    the stepsize by calling takestep.scale() to obtain a certain acceptance
    ration during basing hopping runs.

    Parameters
    ----------
    stepclass : takestep object
        the takestep object which performs the takestep
    acc_ratio : float
        target acceptance ratio
    factor : float
        factor to adjust the stepsize if acc_ratio is too high.
    interval : integer
        adjust the stepsize every interval steps

    note: the keyword frequency is the same as interval.  it exists only for backward compatability
    """

    def __init__(self, stepclass, acc_ratio=0.5, factor=0.9, frequency=None, last_step=None, interval=100,
                 verbose=False):
        self.stepclass = stepclass
        self.accrat = acc_ratio  # target accept ratio
        self.factor = factor
        self.nstepsaccrat = interval
        if frequency is not None:
            print("AdaptiveStepsize: keyword frequency is obsolete, use interval instead")
            self.nstepsaccrat = frequency
        self.last_step = last_step  # stop adjusting after this many steps

        self.naccepted = 0
        self.nsteps = 0
        self.nsteps_tot = 0
        self.verbose = verbose


    def takeStep(self, coords, **kwargs):
        self.nsteps_tot += 1
        self.stepclass.takeStep(coords)

    def updateStep(self, accepted, **kwargs):
        """tell us whether a step was accepted or rejected"""
        self.nsteps += 1
        if self.last_step is not None:
            if self.nsteps_tot > self.last_step:
                return
        if accepted:
            self.naccepted += 1
        if self.nsteps == self.nstepsaccrat:
            self.adjustStep()
            self.nsteps = 0

    def adjustStep(self):
        """adjust the stepsize"""
        rat = float(self.naccepted) / self.nsteps
        if rat > self.accrat:
            self.stepclass.scale(1. / self.factor)
        else:
            self.stepclass.scale(self.factor)

        self.nsteps = 0
        self.naccepted = 0
        if self.verbose:
            print("accrat was ", rat, "new stepsize is ", self.stepclass.stepsize, "f is", self.factor)        

