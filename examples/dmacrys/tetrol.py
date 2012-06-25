import dmagmin_ as GMIN
import pygmin.potentials.gminpotential as gminpot
import numpy as np
import pygmin.basinhopping as bh
from pygmin.optimize import quench
from pygmin.takestep import generic
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
    
def quenchCrystal(coords, pot, **kwargs):
    coords, E, rms, calls = quench.lbfgs_py(coords, pot, **kwargs)
    while(GMIN.reduceCell(coords)):
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
step = DMACRYSTakestep()
step = TestStep()
#0.0
#0.0139495670827
#0.261972630289
#0.356177133898
#0.542180299636
import pickle


if True:
    save = pickle.load(open("storage.1"))
    Eref = save.data[0].E
    GMIN.writeCIF("test.cif", save.data[0].coords)
    for m in save.data:
    #    #print m.coords
        coords, E, tmp, tmp2 = quench.lbfgs_py(m.coords, pot.getEnergyGradient, maxErise=2e-2, M=100, tol=1e-4)
        #print m.E, E
        m.coords = coords
        print m.E - Eref, E - save.data[0].E
        m.E = E
    exit()

save = savenlowest.SaveN(nsave=50, accuracy=1e-4)

opt = bh.BasinHopping(coords, pot, takeStep=step, quenchRoutine=quenchCrystal, temperature=0.4, storage=save)
addColdFusionCheck(opt)
pickle.dump(step, open("bh.dump", "w"))

opt.quenchParameters["tol"]=1e-4
opt.quenchParameters["maxErise"]=2e-2
opt.quenchParameters["maxstep"]=0.1
opt.quenchParameters["M"]=100

for i in xrange(1,50):
    opt.run(100)
    import pickle
    pickle.dump(save, open("storage."+str(i), "w"))

#xa = np.zeros(3*GMIN.getNAtoms())
#GMIN.toAtomistic(xa, coords)
#import pygmin.utils.pymolwrapper as pym
#pym.start()
#pym.draw_spheres(save.da, "A", 1)
