# benchmark all interface
from pele.potentials import LJ

import _pele
import numpy as np
import time
import sys
import _lj
import _lbfgs
from pele.optimize import mylbfgs
import _lj_cython
import _pythonpotential

N=int(sys.argv[2])
natoms=int(sys.argv[1])

print "benchmarking lennard jones potential, %d atoms, %d calls" % (natoms, N)
pot_old = LJ()
pot = _lj.LJ()


t0 = time.time()
for i in xrange(N):
    x = 1.*(np.random.random(3*natoms) - 0.5)
    clbfgs = _lbfgs.LBFGS_CPP(pot, x, tol=1e-4)
    ret = clbfgs.run()

t1 = time.time()
for i in xrange(N):
    x = 1.*(np.random.random(3*natoms) - 0.5)
    ret = mylbfgs(x, pot_old, tol=1e-4)

t2 = time.time()

for i in xrange(N):
    x = 1.*(np.random.random(3*natoms) - 0.5)
    clbfgs = _lbfgs.LBFGS_CPP(_lj_cython.LJ_cython(), x, tol=1e-4)
    ret = clbfgs.run()

t3 = time.time()

for i in xrange(N):
    x = 1.*(np.random.random(3*natoms) - 0.5)
    clbfgs = _lbfgs.LBFGS_CPP(pot_old, x, tol=1e-4)
    ret = clbfgs.run()

t4 = time.time()



print "time for mylbfgs  ", t2-t1
print "time for cpp lbfgs", t1-t0, "speedup",  (t2-t1)/(t1-t0)
print "time for cpp lbfgs with fortran lj", t3-t2, "speedup",  (t2-t1)/(t3-t2)
print "time for cpp lbfgs with old lj", t4-t3, "speedup",  (t2-t1)/(t4-t3)


print ""
print "testing that the potentials all work and give the same result"
potentials = dict(
              lj_old=LJ(),
              lj_cython=_lj_cython.LJ_cython(),
              lj_cpp=_lj.LJ(),
              )
x = np.random.uniform(-1,1,[3*natoms]) * natoms**(1./3)
for name, pot in potentials.iteritems():
    e, g = pot.getEnergyGradient(x)
    label = "%15s" % name
    print label, "getEnergy", pot.getEnergy(x), "getEnergyGradient E", e, "rms", np.linalg.norm(x) / np.sqrt(x.size)

