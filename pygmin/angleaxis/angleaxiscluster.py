import numpy as np
from pygmin.NEB import interpolate_linear
from pygmin.utils.rbtools import CoordsAdapter
from aainterpolate import interpolate_angleaxis

class AngleAxisCluster(CoordsAdapter):
    def __init__(self, nrigid=0, natoms = 0):
        CoordsAdapter.__init__(nrigd=nrigid, natoms=natoms)
    
    def interpolate(self, initial, final, t):
        cinitial = self.copy(initial)
        cfinal = self.copy(final)
        cnew = self.copy(interpolate_linear(self.coords, final, t))
        cnew.rotRigid[:] = interpolate_angleaxis(cinitial.rotRigid, cfinal.rotRigid, t)
        return cnew.coords
    
    def copy(self, coords=None):
        if(coords is None):
            coords = self.coords
        return AngleAxisCluster(nrigid=self.nrigid, natoms=self.natoms, coords=coords)
    
    def distance(self, x1, x2):
        c1 = self.copy(x1)
        c2 = self.copy(x2)
        distance = 0.
        if(not self.posRigid is None):
            distance += np.linalg.norm(self.posRigid)
        if(not self.posAtoms):
            distance += np.linalg.norm(self.posAtoms)
        if(not self.rotRigid):
            distance += np.linalg.norm()