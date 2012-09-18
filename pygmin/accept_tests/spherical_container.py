import numpy as np

class SphericalContainer(object):
    """
    a class to make sure the cluster doesn't leave a sphercial region of
    given radius
    """
    def __init__(self, radius):
        self.radius2 = radius**2
        self.count = 0
        self.nrejected = 0
    
    def accept(self, coords):
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
        return self.accept(coordsnew)
    
    def __call__(self, enew, coordsnew, **kwargs):
        return self.accept(coordsnew)
