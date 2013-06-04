import random
import pickle
import numpy as np
from pele.potentials.gminpotential import GMINPotential
from pele.optimize.quench import mylbfgs, fire
import oxdnagmin_ as GMIN
import time
from pele.systems import oxdna

from optparse import OptionParser
parser = OptionParser()
parser.add_option("--tol", dest="tol", type=float, default=1e-6,
                  help="tolerance for quench")
parser.add_option("-M", dest="M", type=int, default=5,
                  help="M value for mylbfgs")
parser.add_option("--maxErise", dest="maxErise", type=float, default=1e-4,
                  help="maximum energy increase for lbfgs")
parser.add_option("--maxstep", dest="maxstep", type=float, default=0.1,
                  help="maximum stepsize for minimizer")
parser.add_option("-o", dest="outfile", default="results.txt",
                  help="outfile for results")

(options, args) = parser.parse_args()

configurations = pickle.load(open("quench_benchmark.dat"))

GMIN.initialize()
pot = GMINPotential(GMIN)

#fl = open("quenched.xyz", "w")
res = open("results.txt", "w")

t0 = time.time()
print "# energy nsteps rms" 
res.write("# energy nsteps rms\n")
for coords in configurations[0:5]:
    ret = mylbfgs(coords, pot.getEnergyGradient, tol=options.tol, M=options.M, maxstep=options.maxstep, maxErise=options.maxErise)
    energy = ret[1]
    fcalls = ret[3]
    #oxdna.export_xyz(fl, ret[0])
    print energy, fcalls, ret[2]
    res.write("%f %d %g\n"%(energy, fcalls, ret[2]))
              
print "total runtime: ", time.time()-t0
#fl.close()
res.close()