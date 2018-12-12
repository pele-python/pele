from __future__ import print_function
import numpy as np
import pele.exceptions as exc
from . import _spherical_container as fmodule

__all__ = ["SphericalContainer"]


class SphericalContainer(object):
    """
    Reject a structure if any atoms are outside a spherical region

    This test is necessary to avoid evaporation in clusters
    
    a class to make sure the cluster doesn't leave a spherical region of
    given radius.  The center of the spherical region is at the center of mass.
    
    Parameters
    ----------
    radius : float
    """

    def __init__(self, radius, nocenter=False, verbose=False):
        if radius < 0:
            raise exc.SignError
        try:
            self.radius = float(radius)
            self.radius2 = float(radius) ** 2
        except:
            raise TypeError("could not convert {} to float".format(type(radius)))
        self.count = 0
        self.nrejected = 0
        self.nocenter = nocenter
        self.verbose = verbose


    def accept(self, coords):
        """ perform the test"""
        if self.nocenter: return self.accept_fortran(coords)
        self.count += 1
        # get center of mass
        natoms = len(coords) // 3
        coords = np.reshape(coords, [natoms, 3])
        if self.nocenter:
            com = np.zeros(3)
        else:
            com = np.sum(coords, 0) / natoms
        # print np.max(np.sqrt(((coords-com[np.newaxis,:] )**2).sum(1)))
        # print np.max(np.sqrt(((coords)**2).sum(1)))
        reject = (((coords - com[np.newaxis, :] ) ** 2).sum(1) >= self.radius2).any()
        if reject and self.verbose:
            self.nrejected += 1
            print("radius> rejecting", self.nrejected, "out of", self.count)
        return not reject

    def accept_fortran(self, coords):
        return fmodule.check_sphereical_container(coords, self.radius)

    def acceptWrapper(self, eold, enew, coordsold, coordsnew):
        """wrapper for accept"""
        return self.accept(coordsnew)

    def __call__(self, enew, coordsnew, **kwargs):
        """wrapper for accept"""
        return self.accept(coordsnew)

