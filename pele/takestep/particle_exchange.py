from __future__ import print_function
import numpy as np
import random

from pele.takestep.generic import Takestep

__all__ = ["ParticleExchange"]


class ParticleExchange(Takestep):
    """Implement a takestep move which swaps two un-like atoms
    
    Choose a random atom from group A and a random atom from group B
    and exchange their xyz coordinates.
    
    Parameters
    ----------
    Alist, Blist : list of integers
        the indices of the atoms in each of the groups
    verbose : bool
        print debugging info
    """

    def __init__(self, Alist, Blist, verbose=False):
        self.Alist = np.array(Alist)
        self.Blist = np.array(Blist)
        self.verbose = verbose

        self.naccept = 0
        self.ntry = 0

    def takeStep(self, coords, **kwargs):
        iA = random.choice(self.Alist)
        iB = random.choice(self.Blist)
        if self.verbose:
            print("exchange atoms", iA, iB, "accepted", self.naccept, "out of", self.ntry)

        coords = coords.reshape(-1, 3)
        temp = coords[iA, :].copy()
        coords[iA, :] = coords[iB, :]
        coords[iB, :] = temp
        self.ntry += 1
        return coords

    def updateStep(self, accepted, **kwargs):
        """feedback from basin hopping if last step was accepted"""
        if accepted:
            self.naccept += 1
        

