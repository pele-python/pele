import dmagmin_ as GMIN
import pygmin.potentials.gminpotential as gminpot
import pygmin.basinhopping as bh
from pygmin.takestep import generic,group
from pygmin.utils.rbtools import *
from pygmin.potentials.coldfusioncheck import addColdFusionCheck
from pygmin.storage import savenlowest
from pygmin.utils import crystals
import pickle
     
class DMACRYSTakestep(generic.TakestepInterface):
    def takeStep(self, coords, **kwargs):
        GMIN.takestep(coords)

class MyStep(generic.TakestepInterface):
    def takeStep(self, coords, **kwargs):
        from pygmin.takestep import buildingblocks as bb
        ca = CoordsAdapter(nrigid=GMIN.getNRigidBody(), nlattice=6, coords=coords)
        bb.rotate(1.6, ca.rotRigid)
        ca.lattice*=1.2        

def compareMinima(min1, min2):
    from pygmin.utils import rbtools
    ca1 = rbtools.CoordsAdapter(nrigid=GMIN.getNRigidBody(), nlattice=6, coords=min1.coords)
    ca2 = rbtools.CoordsAdapter(nrigid=GMIN.getNRigidBody(), nlattice=6, coords=min2.coords)
    match = crystals.compareStructures(ca1, ca2)
    if not match:
        print "Found minimum with similar energy but different structure"
    return match
    
GMIN.initialize()
pot = gminpot.GMINPotental(GMIN)
crystals.GMIN = GMIN

coords = pot.getCoords()
print "initial energy ", pot.getEnergy(coords)

step1 = MyStep()
reseed = crystals.GenRandomCrystal(CoordsAdapter(nrigid=GMIN.getNRigidBody(), nlattice=6, coords=coords))
step = group.Reseeding(step1, reseed, maxnoimprove=20)

#save = pickle.load(open("storage"))
#save.nsave = 1000
save = savenlowest.SaveN(nsave=100, accuracy=1e-3)
save.compareMinima = compareMinima

opt = bh.BasinHopping(coords, pot, takeStep=step, quenchRoutine=crystals.quenchCrystal, temperature=0.4, storage=save)
addColdFusionCheck(opt)

opt.quenchParameters["tol"]=1e-4
opt.quenchParameters["nsteps"]=1000
opt.quenchParameters["maxErise"]=2e-2
opt.quenchParameters["maxstep"]=0.1
opt.quenchParameters["M"]=100

for i in xrange(1,50):
    opt.run(100)
    pickle.dump(save, open("storage."+str(i), "w"))
    pickle.dump(save, open("storage", "w"))
    i=0
    for m in save.data:
        i+=1
        GMIN.writeCIF("cif/lowest%03d.cif"%(i), m.coords)
