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
nsteps_equil = 10000
nsteps_tot   = 10000

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
    event_after_step=[hist]
    
    radiustest = SphericalContainer(radius)
    accept_tests = [radiustest]
    
    mc = MonteCarlo(coords, potential, takeStep=takestep, temperature=T, \
                    outstream=ostream, event_after_step = event_after_step, \
                    acceptTests = accept_tests)
    mc.printfrq = 1
    replicas.append(mc)


#is it possible to pickle a mc object?
#cp = copy.deepcopy(replicas[0])
#import pickle
#with open("mc.pickle", "w") as fout:
    #pickle.dump(takesteplist[0], fout)



#attach an event to print xyz coords
from printing.print_atoms_xyz import PrintEvent
printxyzlist = []
for n, rep in enumerate(replicas):
    outf = "dumpstruct.%d.xyz" % (n+1) 
    printxyz = PrintEvent(outf, frq=500)
    printxyzlist.append( printxyz)
    rep.addEventAfterStep(printxyz)
    



#do some equilibration steps to find a good step size
for takestep in takesteplist:
    takestep.useAdaptiveStep(last_adaptive_step = nsteps_equil)


ptmc = PTMC(replicas)
ptmc.exchange_frq = 100
ptmc.run(nsteps_tot)



#do production run
#fix the step sizes
#for takestep in takesteplist:
#    takestep.useFixedStep()
#ptmc.run(30000)



print "final energies"
for rep in ptmc.replicas:
    print rep.temperature, rep.markovE
for rep in ptmc.replicas_par:
    print rep.mcsys.markovE
for k in range(nreplicas):
    e,T = ptmc.getRepEnergyT(k)
    print T, e

print "histograms"
for i,hist in enumerate(histograms):
    fname = "hist." + str(i)
    print fname
    with open(fname, "w") as fout:
        for (e, visits) in hist:
            fout.write( "%g %d\n" % (e, visits) )

ptmc.end() #close the open threads
