# -*- coding: iso-8859-1 -*-
import numpy as np
import potentials.lj as lj
#import potentials.ljcpp as lj
import basinhopping as bh
import take_step.random_displacement as random_displacement

natoms = 12
############################################################
#Example 1: Simple basin hopping using defaults
############################################################

# random initial coordinates
coords=np.random.random(3*natoms)

potential = lj.LJ()#1.0, 1.0, None)

step = random_displacement.takeStep( stepsize=0.3)

opt = bh.BasinHopping(coords, potential, takeStep=step.takeStep)
opt.run(100)

############################################################
#Example 2: reading coords from file
############################################################
coords=np.loadtxt('coords')
coords = coords.reshape(coords.size)

step = random_displacement.takeStep( stepsize=0.3)
opt = bh.BasinHopping(coords, potential, takeStep=step.takeStep)

opt.run(100)


############################################################
#Example 3: Saving the coordinates as an xyz file
############################################################
coords=np.random.random(3*natoms)
import saveit

storage = saveit.saveit()
opt = bh.BasinHopping(coords, potential, step.takeStep, storage=storage.insert)
opt.run(100)

with open("lowest", "w") as fout:
  fout.write( str(natoms) + "\n")
  fout.write( str(storage.lowestE) + "\n")
  atoms=storage.lowestcoords.reshape(natoms, 3)  
  for a in atoms: 
      fout.write( "LA "+ str(a[0])+ " " +  str(a[1]) + " " + str(a[2]) + "\n" ) 

############################################################
#Example 4: adaptive step size
############################################################

import take_step.adaptive_step as adaptive_step

manstep = adaptive_step.manageStepSize (0.3, 0.5, 100)
takeStep = random_displacement.takeStep( getStep = manstep.getStepSize )
opt = bh.BasinHopping(coords, potential, takeStep=takeStep.takeStep)
opt.run(100)

