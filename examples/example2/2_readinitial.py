# -*- coding: iso-8859-1 -*-
############################################################
#Example 2: reading coords from file
############################################################
import numpy as np
import potentials.lj as lj
import basinhopping as bh
from takestep import displace

coords=np.loadtxt('coords')
coords = coords.reshape(coords.size)

potential = lj.LJ()

step = displace.RandomDisplacement( stepsize=0.5)
opt = bh.BasinHopping(coords, potential, takeStep=step)

opt.run(100)