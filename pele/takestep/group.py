"""
Created on Jun 7, 2012

@author: vr274
"""
from __future__ import print_function

from pele.takestep import TakestepInterface

__all__ = ["GroupSteps", "BlockMoves", "Reseeding"]


class GroupSteps(TakestepInterface):
    """group several takestep objects

    This class is a wrapper for takestep objects to group several takestep
    objects into one and perform all takesteps simultaneously.

    Parameters
    ----------
    steptakers : list of takestep objects
        list of takestep objects to combine to a single one
    """

    def __init__(self, steptakers):
        self._steptakers = steptakers

    def takeStep(self, coords, **kwargs):
        for step in self._steptakers:
            step.takeStep(coords, **kwargs)

    def updateStep(self, accepted, **kwargs):
        for step in self._steptakers:
            step.updateStep(accepted, **kwargs)


class BlockMoves(TakestepInterface):
    """block based step taking

    This class is a wrapper for takestep objects to group several takestep
    objects into one and perform takestep moves block wise, changing the takestep
    mechanism after a given amount of steps.

    Takestep objects can be added with addBlock

    ::

        step1 = RandomDisplacement()
        step2 = RandomRotation()

        step = BlockMoves()
        step.addBlock(100, step1)
        step.addBlock(100, step2)

    """

    def __init__(self):
        self._steptakers = []
        self._current = 0
        self._counter = 0

    def addBlock(self, nsteps, takestep):
        """add a takestep object

        Parameters
        ----------
        nsteps: integer
            number of steps to perform this type of step
        takestep: takestep object
            takestep object to call for this block

        """
        self._steptakers.append((nsteps, takestep))

    def takeStep(self, coords, **kwargs):
        if self._counter >= self._steptakers[self._current][0]:
            # move to the next step taker 
            self._current += 1
            self._counter = 0
            if self._current >= len(self._steptakers):
                self._current = 0
        self._steptakers[self._current][1].takeStep(coords, **kwargs)
        self._counter += 1

    def updateStep(self, accepted, **kwargs):
        self._steptakers[self._current][1].updateStep(accepted, **kwargs)


class Reseeding(TakestepInterface):
    """Reseeding if energy did not improve

    This class wraps 2 takestep objects, one for a regular takestep move, and a
    reseeding step which is executed if the energy of the basin hopping run did not
    improve after a given number of steps.

    The energy is considered as improved if it is lower than the lowest energy found
    during the current cycle.

    Parameters
    ----------
    takestep : takestep object
        takestep object to perform regular step taking
    reseed : takestep object
        takestep object to do reseeding, e.g. generate a new random configuration
    maxnoimprove : integer
        do 1 reseeding step if the energy did not improve after so many steps.

    """

    def __init__(self, takestep, reseed, maxnoimprove=100, accuracy=1e-4):
        self.takestep = takestep
        self.reseed = reseed
        self.maxnoimprove = maxnoimprove
        self._noimprove = 0
        self.lowest = None
        self.accuracy = accuracy

    def takeStep(self, coords, **kwargs):
        if self._noimprove >= self.maxnoimprove:
            print("The energy did not improve after " + str(self._noimprove) + \
                  " steps, reseeding")
            self.reseed.takeStep(coords, **kwargs)
            kwargs['driver'].acceptTest.forceAccept()
            self.lowest = None
        else:
            self.takestep.takeStep(coords, **kwargs)

    def updateStep(self, accepted, **kwargs):
        driver = kwargs["driver"]
        if self.lowest is None:
            self.lowest = driver.markovE
        if self._noimprove >= self.maxnoimprove:
            self.reseed.updateStep(accepted, **kwargs)
            self._noimprove = 1
        else:
            self.takestep.updateStep(accepted, **kwargs)
            if driver.markovE + self.accuracy >= self.lowest or accepted == False:
                self._noimprove += 1
            else:
                self.lowest = driver.markovE
                self._noimprove = 1
            

