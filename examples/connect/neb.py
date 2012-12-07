import numpy as np
from pygmin.potentials.lj import LJ
from pygmin.transition_states import NEB, InterpolatedPath
from pygmin.mindist import minPermDistStochastic
import pylab as pl
from pygmin import defaults
from copy import copy
from pygmin.optimize import quench

defaults.NEBquenchParams["nsteps"] = 2000
defaults.NEBquenchParams["iprint"] = 1
defaults.NEBquenchParams["maxErise"] = 0.1
defaults.NEBquenchParams["maxstep"] = 0.001
defaults.NEBquenchRoutine = quench.lbfgs_py
#defaults.NEBquenchParams["M"] = 1
k = 100.
nimages = 20
dneb=True

#set up the potential
pot = LJ()

#import the starting and ending points and quench them, 
coords1 = np.genfromtxt("coords.A").flatten()
coords2 = np.genfromtxt("coords.B").flatten()

dist, coords1, coords2 = minPermDistStochastic(coords1, coords2, 
         niter=100, verbose=False)

print "The distance is:", dist

neb = NEB(InterpolatedPath(coords1, coords2, nimages), pot, k=k, dneb=dneb)
neb.optimize()

energies1 = copy(neb.energies)
distances1 = []
for i in xrange(1,len(neb.energies)):
    x1 = neb.coords[i-1]
    x2 = neb.coords[i]
    distances1.append(np.linalg.norm(x2-x1))
neb.optimize()
energies2 = copy(neb.energies)
distances2 = []
for i in xrange(1,len(neb.energies)):
    x1 = neb.coords[i-1]
    x2 = neb.coords[i]
    distances2.append(np.linalg.norm(x2-x1))
    
fig = pl.figure()
pl.plot(energies1)
pl.plot(energies2)
fig = pl.figure()
pl.plot(distances1)
pl.plot(distances2)
pl.show()