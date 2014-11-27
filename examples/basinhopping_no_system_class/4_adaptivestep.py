# -*- coding: iso-8859-1 -*-
# ###########################################################
# Example 4: adaptive step size
############################################################
import numpy as np
import pele.potentials.lj as lj
import pele.basinhopping as bh
from pele.takestep import displace
from pele.takestep import adaptive

natoms = 12

# random initial coordinates
coords = np.random.random(3 * natoms)
potential = lj.LJ()

takeStep = displace.RandomDisplacement(stepsize=0.3)
tsAdaptive = adaptive.AdaptiveStepsize(takeStep, acc_ratio=0.5, frequency=100)
opt = bh.BasinHopping(coords, potential, takeStep=tsAdaptive)
opt.run(1000)