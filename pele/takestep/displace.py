"""
Created on Jun 6, 2012

@author: vr274
"""
from __future__ import absolute_import

import numpy as np
from .generic import TakestepSlice, TakestepInterface
from pele.utils import rotations

__all__ = ["RandomDisplacement", "UniformDisplacement",
           "RotationalDisplacement", "RandomCluster"]


class RandomDisplacement(TakestepSlice):
    """Random displacement on each individual coordinate

    RandomDisplacement is the most basic step taking routine. It simply
    displaces each coordinate my a random value.

    Parameters
    ----------
    stepsize : float
        magnitude of random displacement

    """

    def __init__(self, stepsize=1.0):
        TakestepSlice.__init__(self, stepsize=stepsize)

    def takeStep(self, coords, **kwargs):
        coords[self.srange] += np.random.uniform(low=-self.stepsize, high=self.stepsize, size=coords[self.srange].shape)


class UniformDisplacement(TakestepSlice):
    """Displace each atom be a uniform random vector

    The routine generates a proper uniform random unitvector to displace
    atoms.
    """

    def takeStep(self, coords, **kwargs):
        c = coords[self.srange]
        for x in c.reshape(c.size / 3, 3):
            x += self.stepsize * rotations.vector_random_uniform_hypersphere(3)


class RotationalDisplacement(TakestepSlice):
    """Random rotation for angle axis vector

    RotationalDisplacement performs a proper random rotation. If the coordinate array contains
    positions and orientations, make sure to specify the correct slice for the angle axis
    coordinates.
    """

    def takeStep(self, coords, **kwargs):
        """
        take a random orientational step
        """
        c = coords[self.srange]
        for x in c.reshape(c.size / 3, 3):
            rotations.takestep_aa(x, self.stepsize)


class RandomCluster(TakestepInterface):
    """Generate a random configuration
    """

    def __init__(self, volume=1.0):
        self.volume = volume

    def takeStep(self, coords, **kwargs):
        coords[:] = np.random.random(coords.shape) * (self.volume ** (1. / 3.))
    

