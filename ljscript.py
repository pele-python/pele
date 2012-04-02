# -*- coding: iso-8859-1 -*-
import numpy as np
import potentials.lj as lj
import mingmin
import saveit
import mingmin as gmin
import adaptive_step
import take_step

potential = lj.LJ(1.0, 1.0, None)

coords=np.loadtxt('coords')
# avoid this is this really needed?
natoms=coords.shape[0]
coords = coords.reshape(coords.size)

# run monte carlo
savelowest = saveit.saveit()
manstep = adaptive_step.manageStepSize (0.3, 0.5, 100)
takeStep = take_step.takeStep( RNG = np.random.rand, getStep = manstep.getStepSize )
# potential, coords, natoms, mcsteps, temperature, stepsize, equi steps, savelowest
gmin.monteCarlo(potential, coords, natoms, 1000, 1.0, savelowest, manstep, takeStep)

# print results
with open("lowest", "w") as fout:
  gmin.printxyz(fout,  savelowest.lowestcoords, natoms, savelowest.lowestE)

