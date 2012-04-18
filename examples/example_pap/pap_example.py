'''
Created on Apr 16, 2012

@author: vr274
'''

'''
Created on 13 Apr 2012

@author: ruehle
'''

import pap
import pygmin
import numpy as np
import basinhopping as bh
import take_step.adaptive_step as adaptive_step
import take_step.random_displacement as ts
import storage.savenlowest
        
# initialize GMIN
pygmin.init()

# create PatchyParticle potential
pot = pap.PatchyParticle()    

# load coords array
tmp = np.loadtxt('coords')    
x = tmp.reshape(tmp.size)

# go to reduced coordinate space
x =pot.toReduced(x)

# optional, start with random configuration
x[0:x.size-3] = np.random.random(x.size-3)
    
# use adaptive step size, 0.3 start, acceptance rate 0.5, adjust every 20
manstep = adaptive_step.manageStepSize (0.1, 0.3, 20)    
step = ts.takeStep( getStep=manstep.getStepSize)

# store the lowest 10 minma
minima = storage.savenlowest.SaveN(nsave=10)

# start a basin hopping run
opt = bh.BasinHopping(x, pot,                      
                      temperature=2.,
                      takeStep=step.takeStep,
                      event_after_step=[manstep.insertStepWrapper],
                      storage = minima.insert
                      )

# do 100 mc steps
opt.run(400)

# save the minima
with open("pylowest", "w") as fout:
  for i in minima.data:
      fout.write( str(x.size) + "\n")
      fout.write( "Energy = " + str(i[0]) + "\n")
      for atom in i[1].reshape(x.size/3, 3):
          fout.write(str(atom[0]) + " " + str(atom[1]) + " " + str(atom[2]) + "\n")
