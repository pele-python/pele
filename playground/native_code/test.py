# benchmark all interface
from pele.potentials import LJ

import _pele
import numpy as np
import time
import sys
import _lj
import _lbfgs
from pele.optimize import mylbfgs
N=int(sys.argv[2])
natoms=int(sys.argv[1])

print "benchmarking lennard jones potential, %d atoms, %d calls", natoms, N
pot_old = LJ()
pot = _lj.LJ()


t0 = time.time()
for i in xrange(100):
    x = 1.*(np.random.random(3*natoms) - 0.5)
    clbfgs = _lbfgs.LBFGS_CPP(pot, x, tol=1e-4)
    ret = clbfgs.run()
#e, g = pot.get_energy_gradient(ret[0])
    # print "C", np.linalg.norm(g)

t1 = time.time()
#x = 1.*(np.random.random(3*natoms)-0.5)
for i in xrange(100):
    x = 1.*(np.random.random(3*natoms) - 0.5)
    ret = mylbfgs(x, pot_old, tol=1e-4)
    # print "PY:", np.linalg.norm(pot_old.getEnergyGradient(ret[0])[1])

print "time for mylbfgs  ", time.time()-t1
print "time for cpp lbfgs", t1-t0
print "speedup",  ( time.time()-t1)/(t1-t0)

t0 = time.time()
for i in xrange(N):
    e, g = pot_old.getEnergyGradient(x)
time_f2py = time.time() - t0
print "f2py potential",e,"time",time.time() - t0

t0 = time.time()
for i in xrange(N):
    e, g = pot.getEnergyGradient(x)
time_cython = time.time() - t0
print "cython potential",e,"time",time.time() - t0

