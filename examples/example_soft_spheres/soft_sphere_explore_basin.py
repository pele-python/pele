import numpy as np
from potentials.soft_sphere import SoftSphere, putInBox


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



        


#set up functions to pass to basin hopping

#set up the step taking routine
#Normal basin hopping takes each step from the quenched coords.  This modified step taking routine takes a step from the 
#last accepted coords, not from the quenched coords
from take_step.random_displacement_from_old import takeStep
takestep = takeStep(stepsize=.01)
#this take step routine needs to know if the step was accepted, so pass a function to be called after every step.
event_after_step = [takestep.checkAccepted]

#pass a function which rejects the step if the system leaved the inital basin.
from accept_tests.dont_leave_basin import DontLeaveBasin
accepttest = DontLeaveBasin(1e-4)
accept_test_list = [accepttest.acceptReject]

#pass another function which accepts or rejects based on the energy of the unquenched coords
from accept_tests.metropolis import MetropolisNonQuench
temperature = 1.
mettest = MetropolisNonQuench(temperature, pot)
accept_test_list.append( mettest.acceptReject)
#this accept reject routine also needs to know if the step was accepted, so pass a function to be called after every step.
event_after_step.append( mettest.checkAccepted )


#set up quench routine
#from optimize.quench import fire as quench
#from optimize.quench import cg as quench
from optimize.quench import quench as quench #numpy lbfgs routine
#from optimize.quench import fmin as quench
#from optimize.quench import steepest_descent as quench






#set up basin hopping
from basinhopping import BasinHopping
bh = BasinHopping(coords, pot, takestep.takeStep, \
                  event_after_step = event_after_step, \
                  acceptTests = accept_test_list, \
                  nometropolis = True, \
                  quenchRoutine = quench)

#run basin hopping
bh.run(200)

print bh.naccepted, "steps accepted out of", bh.stepnum
