import pydmagmin as GMIN
import potentials.gminpotential as gminpot
import numpy as np
import basinhopping as bh
from optimize import quench

class DMACRYSTakestep(object):
    def takeStep(self, coords):
        GMIN.takestep(coords)
    
    def updateStep(self, accepted):
        pass

GMIN.initialize()   
pot = gminpot.GMINPotental(GMIN)

coords = pot.getCoords()

step = DMACRYSTakestep()

opt = bh.BasinHopping(coords, pot, takeStep=step, quenchRoutine=quench.lbfgs_py)
opt.run(100)