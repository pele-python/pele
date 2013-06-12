# benchmark all interface
from pele.potentials import LJ

import _pele
import numpy as np
import time
import sys
import _lj
import _lbfgs
from pele.optimize import mylbfgs
N=100
natoms=[10, 13, 20, 30, 31, 38, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150]

print "benchmarking lennard jones potential, %d atoms, %d calls", natoms, N
pot_old = LJ()
pot = _lj.LJ()

clbfgs = _lbfgs.LBFGS(pot)

res = open("results.txt", "w")

for na in natoms:
    t0 = time.time()
    for i in xrange(N):
        x = np.random.random(3*na) - 0.5
        ret = clbfgs.run(x)

    t1 = time.time()
    for i in xrange(N):
        x = np.random.random(3*na) - 0.5
        ret = mylbfgs(x, pot_old, tol=1e-5)
    res.write("%d %f %f\n"%(na, t1-t0, time.time()-t1))
    print "%d %f %f\n"%(na, t1-t0, time.time()-t1)

res.close()

