import numpy as np
from pele.potentials.soft_sphere import SoftSphere, putInBox


natoms = 120
rho = 1.6
boxl = 1.
meandiam = boxl / (float(natoms)/rho)**(1./3)
print "mean diameter", meandiam 

#set up potential
diams = np.array([meandiam for i in range(natoms)]) #make them all the same
pot = SoftSphere(diams = diams)


#initial coordinates
coords = np.random.uniform(-1,1,[natoms*3]) * (natoms)**(1./3)
E = pot.getEnergy(coords)
print "initial energy", E 


#set up quench routine
#from optimize.quench import fire as quench
#from optimize.quench import cg as quench
from pele.optimize import lbfgs_scipy as quench #numpy lbfgs routine
#from optimize.quench import fmin as quench
#from optimize.quench import steepest_descent as quench



#start from quenched coordinates
ret = quench(coords, pot.getEnergyGradient)
coords = ret[0]
        


#set up functions to pass to basin hopping

#set up the step taking routine
#Normal basin hopping takes each step from the quenched coords.  This modified step taking routine takes a step from the 
#last accepted coords, not from the quenched coords
from pele.take_step.random_displacement import takeStep
takestep = takeStep(stepsize=.01)


#pass a function which rejects the step if the system leaved the inital basin.
import do_quenching
dostuff = do_quenching.DoQuenching(pot, coords, quench=quench)
accept_test_list = [dostuff.acceptReject]





#set up basin hopping
from pele.mc import MonteCarlo
temperature = 1.0
event_after_step = []
mc = MonteCarlo(coords, pot, takestep, \
                  event_after_step = event_after_step, \
                  acceptTests = accept_test_list, temperature = temperature)

#run basin hopping
mc.run(200)

print mc.naccepted, "steps accepted out of", mc.stepnum
print "quench: ", dostuff.nrejected, "steps rejected out of", dostuff.ntot

