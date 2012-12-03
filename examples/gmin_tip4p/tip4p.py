from math import sin, cos, pi
import numpy as np
from pygmin.potentials import GMINPotential
from pygmin.angleaxis import RigidFragment, RBSystem
import gmin_ as GMIN
from pygmin.utils import rotations
from pygmin.optimize import quench
from pygmin.takestep import buildingblocks
from pygmin.transition_states import NEB, InterpolatedPath
from pygmin import defaults

# initialize GMIN and potential
GMIN.initialize()
pot = GMINPotential(GMIN)
coords = pot.getCoords()

nrigid = coords.size / 6
print "I have %d water molecules in the system"%nrigid


# define a water molecule
water = RigidFragment()
rho   = 0.9572
theta = 104.52/180.0*pi      
water.add_atom("O", np.array([0., 0., 0.]), 16.)
water.add_atom("H", rho*np.array([0.0, sin(0.5*theta), cos(0.5*theta)]), 1.)
water.add_atom("H", rho*np.array([0.0, -sin(0.5*theta), cos(0.5*theta)]), 1.)
water.finalize_setup()
water.S = np.identity(3, dtype="float64")
# define the whole water system
system = RBSystem()
system.add_sites([water for i in xrange(nrigid)])

# this is an easy access wrapper for coordinates array
ca = system.coords_adapter(coords)

# try rotation symmetry
E = pot.getEnergy(coords)
ca.rotRigid[1] = rotations.rotate_aa(np.array([0., 0., pi]), ca.rotRigid[1])
print "Energy before and after rotation molecule 2 by Pi around its symmetry axis", E, pot.getEnergy(coords)

coords1 = coords.copy()
buildingblocks.rotate(3.0, ca.rotRigid)
ret = quench.mylbfgs(coords, pot.getEnergyGradient, iprint=0)
print ret[1]
coords2 = ret[0]

#np.savetxt("coords1.txt", coords1)
#np.savetxt("coords2.txt", coords2)
coords1 = np.loadtxt("coords1.txt")
coords2 = np.loadtxt("coords2.txt")

defaults.NEBquenchParams["nsteps"] = 500
defaults.NEBquenchParams["iprint"] = 1
defaults.NEBquenchParams["maxErise"] = 0.01
defaults.NEBquenchParams["maxstep"] = 0.05
k = 100.
path=[x for x in InterpolatedPath(coords1, coords2, 30, interpolator=system.interpolate)]
system.align_path(path)

neb = NEB(path, pot, k=k, dneb=True, with_springenergy = True)
neb.optimize()

neb2 = NEB(path, pot, distance=system.neb_distance, k=k/20., dneb=True, with_springenergy = True)
neb2.optimize()
#print path

import pickle
pickle.dump(neb.coords, open("path_cart.dat", "w"))
pickle.dump(neb2.coords, open("path_aa.dat", "w"))

import pylab as pl
fig = pl.figure()
pl.plot(neb.energies, label="Standard NEB")
pl.plot(neb2.energies, label="AngleAxis NEB")
pl.xlabel("neb image")
pl.ylabel("energy")
pl.legend(loc='best')

dist1_cart = []
dist1_aa = []
for i in xrange(0,neb.nimages-1):
    dist1_cart.append(neb.distance(neb.coords[i], neb.coords[i+1])[0])
    dist1_aa.append(system.distance_squared(neb.coords[i], neb.coords[i+1]))

dist2_cart = []
dist2_aa = []
for i in xrange(0,neb2.nimages-1):
    dist2_cart.append(np.linalg.norm(neb2.coords[i] - neb2.coords[i+1])**2)
    dist2_aa.append(system.distance_squared(neb2.coords[i], neb2.coords[i+1]))

dist0_cart = []
dist0_aa = []
for i in xrange(0,neb2.nimages-1):
    dist0_cart.append(np.linalg.norm(path[i] - path[i+1])**2)
    dist0_aa.append(system.distance_squared(path[i], path[i+1]))

fig2 = pl.figure()

pl.plot(dist0_cart, label="interpoaltion - cartesian distances")
pl.plot(dist1_cart, label="Standard NEB - cartesian distances")
pl.plot(dist2_cart, "--", label="AngleAxis NEB - cartesian distances")
pl.xlabel("neb image")
pl.ylabel("energy")
pl.legend(loc='best')
fig3 = pl.figure()

pl.plot(dist0_aa, label="interpoaltion - aa distances")
pl.plot(dist1_aa, label="Standard NEB - aa distances")
pl.plot(dist2_aa, "--", label="AngleAxis NEB - aa distances")
pl.xlabel("neb image")
pl.ylabel("energy")
pl.legend(loc='best')

pl.show()


