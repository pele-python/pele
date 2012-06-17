import dmagmin_ as GMIN
import pygmin.potentials.gminpotential as gminpot
import numpy as np
import pygmin.basinhopping as bh
from pygmin.optimize import quench
from pygmin.takestep import generic
from pygmin.utils.rbtools import *
from pygmin.potentials.coldfusioncheck import addColdFusionCheck

class DMACRYSTakestep(generic.TakestepInterface):
    def takeStep(self, coords, **kwargs):
        GMIN.takestep(coords)

class TestStep(generic.TakestepInterface):
    def takeStep(self, coords, **kwargs):
        from pygmin.takestep import buildingblocks as bb
        print coords.size
        ca = CoordsAdapter(nrigid=2, nlattice=6, coords=coords)
        bb.rotate(1.6, ca.rotRigid)
        ca.lattice*=1.1
    
GMIN.initialize()   
pot = gminpot.GMINPotental(GMIN)

coords = pot.getCoords()

#print quench.lbfgs_py(coords, pot.getEnergyGradient)

#exit()
step = DMACRYSTakestep()
step = TestStep()

opt = bh.BasinHopping(coords, pot, takeStep=step, quenchRoutine=quench.lbfgs_py, temperature=0.1)
addColdFusionCheck(opt)
opt.quenchParameters["tol"]=1e-3
opt.quenchParameters["maxErise"]=2e-2
opt.quenchParameters["maxstep"]=0.1
opt.quenchParameters["M"]=100

opt.run(100)