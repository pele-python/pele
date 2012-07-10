import dmagmin_ as GMIN
import pygmin.potentials.gminpotential as gminpot
import numpy as np
import pygmin.basinhopping as bh
from pygmin.optimize import quench
from pygmin.takestep import generic,group
from pygmin.utils.rbtools import *
from pygmin.potentials.coldfusioncheck import addColdFusionCheck
from pygmin.storage import savenlowest
     
class DMACRYSTakestep(generic.TakestepInterface):
    def takeStep(self, coords, **kwargs):
        GMIN.takestep(coords)

class TestStep(generic.TakestepInterface):
    def takeStep(self, coords, **kwargs):
        from pygmin.takestep import buildingblocks as bb
        print coords.size
        ca = CoordsAdapter(nrigid=2, nlattice=6, coords=coords)
        bb.rotate(1.6, ca.rotRigid)
        ca.lattice*=1.2
        
class GenRandomCrystal(generic.TakestepInterface):
    def __init__(self, volume=None, shear=2., expand=2.0):
        self.volume = volume
        self.shear = shear
        self.expand = expand
        
    def takeStep(self, coords, **kwargs):
        from pygmin import rotations
        from pygmin.utils import lattice
        ca = CoordsAdapter(nrigid=2, nlattice=6, coords=coords)
        
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

        
def quenchCrystal(coords, pot, **kwargs):
    coords, E, rms, calls = quench.lbfgs_py(coords, pot, **kwargs)
    #while(GMIN.reduceCell(coords)):
    GMIN.reduceCell(coords)
    print "Reduced cell, redo minimization"
    coords, E, rms, callsn = quench.lbfgs_py(coords, pot, **kwargs)
    calls+=callsn
    return coords, E, rms, calls

GMIN.initialize()   
pot = gminpot.GMINPotental(GMIN)
print "get coords"
coords = pot.getCoords()
print "initial energy ", pot.getEnergy(coords)
#print GMIN.getNRigidBody()
#print quench.lbfgs_py(coords, pot.getEnergyGradient)
#print GMIN.reduceCell(coords)

# 1.43174012081
#exit()
step1 = TestStep()
reseed = GenRandomCrystal()
step = group.Reseeding(step1, reseed, maxnoimprove=20)

#0.0
#0.0139495670827
#0.261972630289
#0.356177133898
#0.542180299636
import pickle

print "hallo"

save = pickle.load(open("storage"))
save.nsave = 1000
#save = savenlowest.SaveN(nsave=50, accuracy=1e-4)

if False:
    Eref = save.data[0].E
    GMIN.writeCIF("test.cif", save.data[0].coords)
    print "blas"
    i=1
    print save
    for m in save.data[0:2]:
        #        m = save.data[0]
        #    #print m.coords
        coords, E, tmp, tmp2 = quenchCrystal(m.coords, pot.getEnergyGradient, maxErise=2e-2, M=100, tol=1e-10)
        #print m.E, E
        #m.coords = coords
        #print m.E - Eref, E - save.data[0].E
        #print m.E - Eref
        print m.E
        #m.E = E
        print "Exporting"
        #print GMIN.reduceCell(m.coords)
        #print "lowest%3d.cif"%(i)
        GMIN.writeCIF("cif/lowest%03d.cif"%(i), coords)
        print "do structut", i
        i+=1

    exit()

opt = bh.BasinHopping(coords, pot, takeStep=step, quenchRoutine=quenchCrystal, temperature=0.4, storage=save)
addColdFusionCheck(opt)
pickle.dump(step, open("bh.dump", "w"))

opt.quenchParameters["tol"]=1e-4
opt.quenchParameters["nsteps"]=1000
opt.quenchParameters["maxErise"]=2e-2
opt.quenchParameters["maxstep"]=0.1
opt.quenchParameters["M"]=100

for i in xrange(1,50):
    opt.run(100)
    import pickle
    pickle.dump(save, open("storage."+str(i), "w"))
    pickle.dump(save, open("storage", "w"))
    i=0
    for m in save.data:
        i+=1
        GMIN.writeCIF("cif/lowest%03d.cif"%(i), m.coords)

#xa = np.zeros(3*GMIN.getNAtoms())
#GMIN.toAtomistic(xa, coords)
#import pygmin.utils.pymolwrapper as pym
#pym.start()
#pym.draw_spheres(save.da, "A", 1)
