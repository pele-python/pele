import numpy as np
import potentials.lj as lj
#import potentials.ljcpp as lj
from mc import MonteCarlo 
import take_step.random_displacement as random_displacement
from ptmc import PTMC, getTemps
import copy
from tools.histogram import EnergyHistogram
from optimize.quench import quench

natoms = 19

coords=np.random.random(3*natoms)
#quench the coords so we start from a reasonable location
mypot = lj.LJ()
ret = quench(coords, mypot.getEnergyGradient)
coords = ret[0]

nreplicas = 4
Tmin = 0.2
Tmax = 0.4

Tlist = getTemps(Tmin, Tmax, nreplicas)
replicas = []
ostreams = []
histograms = []
for i in range(nreplicas):
    T = Tlist[i]
    potential = lj.LJ()
    takestep = random_displacement.takeStep( stepsize=0.1)
    file = "mcout." + str(i+1)
    ostream = open(file, "w")
    hist = EnergyHistogram( -60., 60., 100)
    histograms.append(hist)
    event_after_step=[hist.insertWrapper]
    mc = MonteCarlo(coords, potential, takeStep=takestep, temperature=T, \
                    outstream=ostream, event_after_step = event_after_step)
    replicas.append(mc)


ptmc = PTMC(replicas)
ptmc.exchange_frq = 10000000
ptmc.run(100000)

print "final energies"
for rep in replicas:
    print rep.temperature, rep.markovE

print "histograms"
for i,hist in enumerate(histograms):
    fname = "hist." + str(i)
    print fname
    with open(fname, "w") as fout:
        for (e, visits) in hist:
            fout.write( "%g %d\n" % (e, visits) )