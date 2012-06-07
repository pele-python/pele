# -*- coding: iso-8859-1 -*-
import numpy as np
import potentials.lj as lj
#import potentials.ljcpp as lj
import basinhopping as bh
from takestep import displace

natoms = 12
############################################################
#Example 1: Simple basin hopping using defaults
############################################################

# random initial coordinates
coords=np.random.random(3*natoms)

potential = lj.LJ()#1.0, 1.0, None)

step = displace.RandomDisplacement(stepsize=0.3)

opt = bh.BasinHopping(coords, potential, takeStep=step)
opt.run(100)

############################################################
#Example 2: reading coords from file
############################################################

coords=np.loadtxt('coords')
coords = coords.reshape(coords.size)

step = displace.RandomDisplacement( stepsize=0.3)
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
  fout.write( str(storage.data[0].E) + "\n")
  atoms=storage.data[0].coords.reshape(natoms, 3)  
  for a in atoms: 
      fout.write( "LA "+ str(a[0])+ " " +  str(a[1]) + " " + str(a[2]) + "\n" ) 

############################################################
#Example 4: adaptive step size
############################################################

from takestep import adaptive

takeStep = displace.RandomDisplacement( stepsize=0.3 )
tsAdaptive = adaptive.AdaptiveStepsize(takeStep, acc_ratio = 0.5, frequency = 100)
opt = bh.BasinHopping(coords, potential, takeStep=tsAdaptive)
opt.run(1000)

