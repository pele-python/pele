# -*- coding: iso-8859-1 -*-
import numpy as np
import potentials.lj as lj
import basinhopping as bh
from takestep import displace
import pickle

print "Setting up basin hopping class"

natoms = 12
coords=np.random.random(3*natoms)
potential = lj.LJ()#1.0, 1.0, None)
step = displace.CartesianDisplacement(stepsize=0.3)
opt = bh.BasinHopping(coords, potential, takeStep=step)

print "Do 10 basin hopping steps"
opt.run(10)

#try:
print "Save state to file bh.dump"
output = open("bh.dump", "w")
pickle.dump(opt, output)
del output
#except TypeError as err:
#    print err

print "Continue with 10 steps from file"
infile = open("bh.dump", "r")
opt2 = pickle.load(infile)
opt2.potential = potential
opt2.run(10)