import numpy as np
import pygmin.exceptions as exc

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
    def __init__(self, radius):
        if radius < 0:
            raise exc.SignError
        self.radius2 = float(radius)**2
        self.count = 0
        self.nrejected = 0
    
    def accept(self, coords):
        """ perform the test"""
        self.count += 1
        #get center of mass
        natoms = len(coords)/3
        coords = np.reshape(coords, [natoms,3])
        com = np.sum(coords, 0)/natoms
        reject = ( ((coords-com[np.newaxis,:] )**2).sum(1) >= self.radius2 ).any()
        if reject: 
            self.nrejected += 1
            print "radius> rejecting", self.nrejected, "out of", self.count
        return not reject
    
    def acceptWrapper(self, eold, enew, coordsold, coordsnew):
        """wrapper for accept"""
        return self.accept(coordsnew)
    
    def __call__(self, enew, coordsnew, **kwargs):
        """wrapper for accept"""
        return self.accept(coordsnew)
