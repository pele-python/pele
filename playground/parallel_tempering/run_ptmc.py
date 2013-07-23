import numpy as np
import pele.potentials.lj as lj
#import potentials.ljcpp as lj
from pele.mc import MonteCarlo 
from pele.takestep import RandomDisplacement, AdaptiveStepsize
from ptmc import PTMC, getTemps
import copy
from pele.utils.histogram import EnergyHistogram, PrintHistogram
from pele.optimize import mylbfgs as quench
from pele.accept_tests.spherical_container import SphericalContainer

def runptmc(nsteps_tot = 100000):
    natoms = 31
    nreplicas = 4
    Tmin = 0.2
    Tmax = 0.4
    
    nsteps_equil = 10000
    nsteps_tot   = 100000
    histiprint  =  nsteps_tot / 10
    exchange_frq = 100*nreplicas
    
    coords=np.random.random(3*natoms)
    #quench the coords so we start from a reasonable location
    mypot = lj.LJ()
    ret = quench(coords, mypot)
    coords = ret.coords
    
    
    Tlist = getTemps(Tmin, Tmax, nreplicas)
    replicas = []
    ostreams = []
    histograms = []
    takesteplist = []
    radius = 2.5
    # create all the replicas which will be passed to PTMC
    for i in range(nreplicas):
        T = Tlist[i]
        potential = lj.LJ()
        
        takestep = RandomDisplacement( stepsize=0.01)
        adaptive = AdaptiveStepsize(takestep, last_step = nsteps_equil)
        takesteplist.append( adaptive )
        
        file = "mcout." + str(i+1)
        ostream = open(file, "w")
        
        hist = EnergyHistogram( -134., 10., 1000)
        histograms.append(hist)
        event_after_step=[hist]
        
        radiustest = SphericalContainer(radius)
        accept_tests = [radiustest]
        
        mc = MonteCarlo(coords, potential, takeStep=takestep, temperature=T, \
                        outstream=ostream, event_after_step = event_after_step, \
                        confCheck = accept_tests)
        mc.histogram = hist #for convienence
        mc.printfrq = 1
        replicas.append(mc)
    
    
    #is it possible to pickle a mc object?
    #cp = copy.deepcopy(replicas[0])
    #import pickle
    #with open("mc.pickle", "w") as fout:
        #pickle.dump(takesteplist[0], fout)
    
    
    
    #attach an event to print xyz coords
    from pele.printing.print_atoms_xyz import PrintEvent
    printxyzlist = []
    for n, rep in enumerate(replicas):
        outf = "dumpstruct.%d.xyz" % (n+1) 
        printxyz = PrintEvent(outf, frq=500)
        printxyzlist.append( printxyz)
        rep.addEventAfterStep(printxyz)
        
    
    #attach an event to print histograms
    for n, rep in enumerate(replicas):
        outf = "hist.%d" % (n+1)
        histprint = PrintHistogram(outf, rep.histogram, histiprint)
        rep.addEventAfterStep(histprint)
    
    
    ptmc = PTMC(replicas)
    ptmc.use_independent_exchange = True
    ptmc.exchange_frq = exchange_frq
    ptmc.run(nsteps_tot)
    
    
    
    #do production run
    #fix the step sizes
    #for takestep in takesteplist:
    #    takestep.useFixedStep()
    #ptmc.run(30000)
    
    
    if False: #this doesn't work
        print "final energies"
        for rep in ptmc.replicas:
            print rep.temperature, rep.markovE
        for rep in ptmc.replicas_par:
            print rep.mcsys.markovE
        for k in range(nreplicas):
            e,T = ptmc.getRepEnergyT(k)
            print T, e
    
    if False: #this doesn't work
        print "histograms"
        for i,hist in enumerate(histograms):
            fname = "hist." + str(i)
            print fname
            with open(fname, "w") as fout:
                for (e, visits) in hist:
                    fout.write( "%g %d\n" % (e, visits) )
    
    ptmc.end() #close the open threads
    
def getReplicaPath(fname = "exchanges", nreps = 4):
    paths = [ [i] for i in range(nreps)]
    positions = np.array(range(nreps))
    newpositions = np.array(range(nreps))
    oldpositions = np.array(range(nreps))
    with open(fname, "r") as fin:
        for line in fin:
            sline = line.split()
            time = int(sline[0])
            for newposition in range(nreps):
                oldposition = int(sline[3+newposition])
                oldpositions[newposition] = oldposition
                #replica = position2replica[oldposition]
                #newpositions[replica] =
            print oldpositions 
            print positions
            #positions[:] = positions[newpositions]
            positions[:] = positions[oldpositions]
            print positions
            #print ""
            for i, j in enumerate(positions):
                paths[i].append(j)
    if True:
        import matplotlib.pyplot as plt
        nppaths = np.array(paths)
        print np.shape(nppaths)
        for i in range(4):
            plt.subplot(2,2,i+1)
            plt.plot( nppaths[i,:])
        plt.show()
        print paths[0]
            


if __name__ == "__main__":
    if True:
        runptmc(
                nsteps_tot = 1000000
                )
    getReplicaPath()
