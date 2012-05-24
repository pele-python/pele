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

opt = bh.BasinHopping(coords, potential, takeStep=step)
opt.run(100)

############################################################
#Example 2: reading coords from file
############################################################
coords=np.loadtxt('coords')
coords = coords.reshape(coords.size)

step = random_displacement.takeStep( stepsize=0.3)
opt = bh.BasinHopping(coords, potential, takeStep=step)

opt.run(200)


############################################################
#Example 3: Saving the coordinates as an xyz file
############################################################
coords=np.random.random(3*natoms)
from storage.savenlowest import SaveN as saveit

storage = saveit()
opt = bh.BasinHopping(coords, potential, step, storage=storage.insert)
opt.run(200)

with open("lowest", "w") as fout:
  fout.write( str(natoms) + "\n")
  fout.write( str(storage.data[0][0]) + "\n")
  atoms=storage.data[0][1].reshape(natoms, 3)  
  for a in atoms: 
      fout.write( "LA "+ str(a[0])+ " " +  str(a[1]) + " " + str(a[2]) + "\n" ) 

############################################################
#Example 4: adaptive step size
############################################################

takeStep = random_displacement.takeStep( stepsize=0.3 )
takeStep.useAdaptiveStep(acc_ratio = 0.5, freq = 100)
opt = bh.BasinHopping(coords, potential, takeStep=takeStep)
opt.run(200)

