import dmagmin_ as GMIN
from pygmin.utils import rotations
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
    ''' takestep class to generate a random crystal
    
        GenRandomCrystal is a takestep class which generates a random crystal structure.
        It can be either used as a standard takestep routine to perform a random search
        or as a reseeding routine in combination with pygmin.takestep.group.Reseeding
    '''
    
    def __init__(self, coordsadapter, volume=None, shear=2., expand=2.0):
        '''
        :param volume: volume bounds of the generated cell [min, max]. None to use 2* current volume
        :param shear: maximum off diagonal matrix element for lattice matrix
        :param expane: maxumum assymmetry of the cell, 0 means generate a cubic cell
        
        ''' 
        
        self.volume = volume
        self.shear = shear
        self.expand = expand
        self.coordsadapter = coordsadapter
        
    def takeStep(self, coords, **kwargs):
        ''' takeStep routine to generate random cell '''        
        ca = self.coordsadapter        
        ca.updateCoords(coords)
        
        volumeTarget = 2.*lattice.volume(ca.lattice)
        # first choose random positions and rotations
        for i in xrange(GMIN.getNRigidBody()):
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
    ''' Special quench routine for crystals which makes sure that the final structure is a reduced cell '''
    coords, E, rms, calls = quench.lbfgs_py(coords, pot, **kwargs)
    #while(GMIN.reduceCell(coords)):
    if(GMIN.reduceCell(coords)):
        #print "Reduced cell, redo minimization"
        coords, E, rms, callsn = quench.lbfgs_py(coords, pot, **kwargs)
        calls+=callsn
    return coords, E, rms, calls        