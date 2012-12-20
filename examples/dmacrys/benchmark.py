import dmagmin_ as GMIN
import pygmin.potentials.gminpotential as gminpot
import numpy as np
import pygmin.basinhopping as bh
import pygmin.optimize._quench as quench
from pygmin.takestep import generic
from pygmin.utils.rbtools import *
from pygmin.potentials.coldfusioncheck import addColdFusionCheck

GMIN.initialize()   
#pot = PotWrap(GMIN) 
pot = gminpot.GMINPotential(GMIN)

coords = pot.getCoords()
print "initial energy is a", pot.getEnergy(coords)

from pygmin.utils import benchmark

print "quenching to minimum"
coords, E, tmp, tmp=quench.lbfgs_py(coords, pot.getEnergyGradient, tol=1e-3, maxErise=2e-4)
    
print "do small displacement"
coords = coords + np.random.random(coords.shape)*0.01
print "check new minimum (hopefully the same)"
tmp, Emin, tmp, tmp = quench.lbfgs_py(coords, pot.getEnergyGradient, tol=1e-3)

print "setting up benchmark"
bench = benchmark.QuenchBenchmark(pot)
#bench.addMinimizer("lbfgs", quench.quench)
#bench.addMinimizer("mylbfgs", quench.mylbfgs)
bench.addMinimizer("lbfgs_py", quench.lbfgs_py)
#bench.addMinimizer("lbfgs_ase", quench.lbfgs_ase)
#bench.addMinimizer("cg", quench.cg)
bench.addMinimizer("fire", quench.fire)
#bench.addMinimizer("bfgs", quench.bfgs)
#bench.addMinimizer("fmin", quench.fmin)
#bench.addMinimizer("steep", quench.steepest_descent)

print "The reference energy is " + str(Emin)
bench.run(Emin,coords)
bench.plot()
