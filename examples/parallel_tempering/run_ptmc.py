import numpy as np
import potentials.lj as lj
#import potentials.ljcpp as lj
from mc import MonteCarlo 
import take_step.random_displacement as random_displacement
from ptmc import PTMC, getTemps
import copy

natoms = 19

coords=np.random.random(3*natoms)
nreplicas = 4
Tmin = 1.
Tmax = 1.2

Tlist = getTemps(Tmin, Tmax, nreplicas)
replicas = []
ostreams = []
for i in range(nreplicas):
    T = Tlist[i]
    potential = lj.LJ()
    takestep = random_displacement.takeStep( stepsize=0.1)
    file = "mcout." + str(i+1)
    ostream = open(file, "w")
    mc = MonteCarlo(coords, potential, takeStep=takestep, temperature=T, outstream=ostream)
    replicas.append(mc)


ptmc = PTMC(replicas)
ptmc.run(1000)