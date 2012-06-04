import pydmagmin as GMIN
import potentials.gminpotential as gminpot
import numpy as np
import basinhopping as bh
import take_step.random_displacement as random_displacement
from optimize import quench

class DMACRYSTakestep(object):
    def takeStep(self, coords):
        GMIN.takestep(coords)
    
    def updateStep(self, accepted):
        pass

GMIN.initialize()   
pot = gminpot.GMINPotental(GMIN)
coords = pot.getCoords()

print GMIN.getDOF(), GMIN.getNAtoms()
print coords
print pot.getEnergyGradient(coords)

step = DMACRYSTakestep()

opt = bh.BasinHopping(coords, pot, takeStep=step, quenchRoutine=quench.fire)
opt.run(100)