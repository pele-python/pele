# benchmark all interface
from pygmin.potentials import LJ

import _pygmin
import numpy as np
import time
import sys
import _lj

N=int(sys.argv[2])
natoms=int(sys.argv[1])

print "benchmarking lennard jones potential, %d atoms, %d calls", natoms, N
pot_old = LJ(0.5, 1.0)
pot = _lj.LJ(4.*0.5, 4.*0.5 )


x = np.random.random(3*natoms)

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

