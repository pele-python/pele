import dmagmin_ as GMIN
from pygmin import rotations
from pygmin.takestep import generic
from pygmin.optimize import quench
import crystals
import lattice
import numpy as np

# compare 2 minima
def compareMinima(min1, min2):
    from pygmin.utils import rbtools
    ca1 = rbtools.CoordsAdapter(nrigid=GMIN.getNRigidBody(), nlattice=6, coords=min1.coords)
    ca2 = rbtools.CoordsAdapter(nrigid=GMIN.getNRigidBody(), nlattice=6, coords=min2.coords)
    match = crystals.compareStructures(ca1, ca2)
    return match

# generate a random crystal structure
class GenRandomCrystal(generic.TakestepInterface):
    def __init__(self, coordsadapter, volume=None, shear=2., expand=2.0):
        self.volume = volume
        self.shear = shear
        self.expand = expand
        self.coordsadapter = coordsadapter
        
    def takeStep(self, coords, **kwargs):        
        ca = self.coordsadapter        
        ca.updateCoords(coords)
        
        volumeTarget = 2.*lattice.volume(ca.lattice)
        # first choose random positions and rotations
        for i in xrange(2):
            ca.posRigid[i] = np.random.random()
            ca.rotRigid[i] = rotations.random_aa()
         
        # random box
        ca.lattice[[0,3,5]] = 1.0 + self.expand * np.random.random(3)  
        ca.lattice[[1,2,4]] = self.shear * np.random.random(3)
        
        if(self.volume != None):
            volumeTarget = self.volume[0] + (self.volume[1] - self.volume[0]) * np.random.random()
                    
        vol = lattice.volume(ca.lattice)
        ca.lattice[:] = ca.lattice * (volumeTarget / vol)**(1.0/3.0)
        GMIN.reduceCell(coords)
        
# special quencher for crystals
def quenchCrystal(coords, pot, **kwargs):
    coords, E, rms, calls = quench.lbfgs_py(coords, pot, **kwargs)
    #while(GMIN.reduceCell(coords)):
    if(GMIN.reduceCell(coords)):
        #print "Reduced cell, redo minimization"
        coords, E, rms, callsn = quench.lbfgs_py(coords, pot, **kwargs)
        calls+=callsn
    return coords, E, rms, calls        