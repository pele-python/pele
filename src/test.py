from pygmin.potentials import LJ

import _pygmin
import numpy as np
import time
import sys
import _lj
from _lj_cython import LJ_cython

N=int(sys.argv[2])
natoms=int(sys.argv[1])

sigma = 0.5
epsilon = 0.8

print "benchmarking lennard jones potential, %d atoms, %d calls", natoms, N
pot_old = LJ(sig=sigma, eps=epsilon)
pot = _lj.LJ(sigma=sigma, epsilon=epsilon)
pot2 = LJ_cython(sigma=sigma, epsilon=epsilon)
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

grad = x.copy()
t0 = time.time()
for i in xrange(N):
    e = pot.get_energy_gradient_inplace(x, grad)
time_inplace = time.time() - t0
print "cython in place",e,"time",time.time() - t0

t0 = time.time()
_pygmin.call_pot(pot, x, grad, N)
time_cloop = time.time() - t0
print "loop in c time",time.time() - t0

t0 = time.time()
for i in xrange(N):
    e, g = pot2.get_energy_gradient(x)
time_cython2 = time.time() - t0
print "cython LJ potential",e,"time",time.time() - t0

t0 = time.time()
_pygmin.call_pot(pot2, x, grad, N)
time_cloop2 = time.time() - t0
print "cython LJ potentiai cloop",e,"time",time.time() - t0

print 
print "# natoms N f2py cython inplace cloop cython2 cloop2"
print natoms, N, time_f2py, time_cython, time_inplace, time_cloop, time_cython2, time_cloop2
