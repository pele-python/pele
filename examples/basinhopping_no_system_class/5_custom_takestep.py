# -*- coding: iso-8859-1 -*-
# ###########################################################
# Example 5: Adding a custom takestep routine.  This example
# takes 100 monte carlo steps as one basin hopping step
# ###########################################################
import numpy as np
import pele.potentials.lj as lj
import pele.basinhopping as bh
from pele.takestep import displace
from pele.mc import MonteCarlo


class TakeStepMonteCarlo:
    def __init__(self, pot, T=10., nsteps=100, stepsize=0.1):
        self.potential = pot
        self.T = T
        self.nsteps = nsteps

        self.mcstep = displace.RandomDisplacement(stepsize=stepsize)

    def takeStep(self, coords, **kwargs):
        # make a new monte carlo class
        mc = MonteCarlo(coords, self.potential, self.mcstep,
                        temperature=self.T, outstream=None)
        mc.run(self.nsteps)
        coords[:] = mc.coords[:]

    def updateStep(self, acc, **kwargs):
        pass


natoms = 12

# random initial coordinates
coords = np.random.random(3 * natoms)
potential = lj.LJ()

step = TakeStepMonteCarlo(potential)

opt = bh.BasinHopping(coords, potential, takeStep=step)
opt.run(100)

# some visualization
try:
    import pele.utils.pymolwrapper as pym

    pym.start()
    pym.draw_spheres(opt.coords, "A", 1)
except:
    print "Could not draw using pymol, skipping this step"