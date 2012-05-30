import numpy as np
import potentials.lj as lj
#import potentials.ljcpp as lj
from mc import MonteCarlo 
import take_step.random_displacement as random_displacement
from ptmc import PTMC, getTemps
import copy
from tools.histogram import EnergyHistogram
from optimize.quench import quench
from accept_tests.spherical_container import SphericalContainer

natoms = 31

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
takesteplist = []
radius = 2.5
for i in range(nreplicas):
    T = Tlist[i]
    potential = lj.LJ()
    takestep = random_displacement.takeStep( stepsize=0.01)
    takesteplist.append( takestep )
    file = "mcout." + str(i+1)
    ostream = open(file, "w")
    hist = EnergyHistogram( -134., 10., 1000)
    histograms.append(hist)
    event_after_step=[hist.insertWrapper]
    
    radiustest = SphericalContainer(radius)
    accept_tests = [radiustest.acceptWrapper]
    
    mc = MonteCarlo(coords, potential, takeStep=takestep, temperature=T, \
                    outstream=ostream, event_after_step = event_after_step, \
                    acceptTests = accept_tests)
    mc.printfrq = 1000
    replicas.append(mc)


#attach an event to print xyz coords
from printing.print_atoms_xyz import PrintEvent
printxyzlist = []
for n, rep in enumerate(replicas):
    outf = "dumpstruct.%d.xyz" % (n+1) 
    printxyz = PrintEvent(outf, frq=500)
    printxyzlist.append( printxyz)
    rep.addEventAfterStep(printxyz.event)


ptmc = PTMC(replicas)
ptmc.exchange_frq = 100

#do some equilibration steps to find a good step size
for takestep in takesteplist:
    takestep.useAdaptiveStep()
ptmc.run(10000)

#do production run
#fix the step sizes
for takestep in takesteplist:
    takestep.useFixedStep()
ptmc.run(10000000)

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
