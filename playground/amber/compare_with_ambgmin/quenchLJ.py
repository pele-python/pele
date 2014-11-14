############################################################
# Just quench 
############################################################
import numpy as np

import pele.potentials.lj as lj

natoms = 12

# random initial coordinates
coords=np.random.random(3*natoms)
pot = lj.LJ()

print pot.getEnergy(coords) 

a,b = pot.getEnergyGradient(coords) 
print type(a) 

from pele.optimize import lbfgs_scipy, cg  fire

# lbfgs 
ret = lbfgs_scipy( coords, pot, iprint=-1 , tol = 1e-3, nsteps=100) 
     # core dump! 

# cg  
#    ret = cg( coordsVec, pot.getEnergyGradient) 
    # runtime error -- ValueError: The truth value of an array with more than ...
    
# fire   
#ret = fire( coords, pot.getEnergyGradient, tol = 1e-3, nsteps=20000) 
# ValueError: The truth value of an array with more than ...
    # works but after 1000 iterations gives an energy of -90.9378267921 higher than initial energy of -90.9364375726!
    
print "energy ", ret.energy
print "rms gradient", ret.rms
print "number of function calls", ret.nfev

