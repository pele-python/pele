# benchmark all interface
from pygmin.potentials import LJ

import _pygmin
import numpy as np
import time
import sys
import _lj
import _lbfgs
from pygmin import optimize
N=10 # int(sys.argv[2])
natoms=38 #int(sys.argv[1])

print "benchmarking lennard jones potential, %d atoms, %d calls", natoms, N
pot_old = LJ(0.5, 1.0)
pot = _lj.LJ(4.*0.5, 4.*0.5 )

clbfgs = _lbfgs.LBFGS(pot)

t0 = time.time()
for i in xrange(100):
    x = np.random.random(3*natoms)
    ret = clbfgs.run(x)
    print "C", np.linalg.norm(pot.get_energy_gradient(ret[0])[1])

t1 = time.time()
for i in xrange(100):
    x = np.random.random(3*natoms)
    ret = optimize.mylbfgs(x, pot_old.getEnergyGradient, tol=1e-5)
    print "PY:", np.linalg.norm(pot_old.getEnergyGradient(ret[0])[1])

print time.time()-t1, t1-t0

t0 = time.time()
for i in xrange(N):
    e, g = pot_old.getEnergyGradient(x)
time_f2py = time.time() - t0
print "f2py potential",e,"time",time.time() - t0

t0 = time.time()
for i in xrange(N):
    e, g = pot.get_energy_gradient(x)
time_cython = time.time() - t0
print "cython potential",e,"time",time.time() - t0

