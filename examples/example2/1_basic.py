# -*- coding: iso-8859-1 -*-
############################################################
#Example 1: Simple basin hopping using defaults
############################################################
import numpy as np
import potentials.lj as lj
import basinhopping as bh
from takestep import displace

natoms = 12

# random initial coordinates
coords=np.random.random(3*natoms)
potential = lj.LJ()

step = displace.RandomDisplacement(stepsize=0.5)

opt = bh.BasinHopping(coords, potential, takeStep=step)
opt.run(100)

# some visualization
try: 
    import utils.pymolwrapper as pym
    pym.start()
    pym.draw_spheres(opt.coords, "A", 1)
except:
    print "Could not draw using pymol, skipping this step"