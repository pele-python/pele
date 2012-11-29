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

defaults.NEBquenchParams["nsteps"] = 100
defaults.NEBquenchParams["iprint"] = 1


neb = NEB(InterpolatedPath(coords1, coords2, 50), pot)
neb.optimize()

neb2 = NEB(InterpolatedPath(coords1, coords2, 50), pot, distance=system.neb_distance)
neb2.optimize()

import pylab as pl
pl.plot(neb.energies, label="Standard NEB")
pl.plot(neb2.energies, label="AngleAxis NEB")
pl.xlabel("neb image")
pl.ylabel("energy")
pl.legend()
pl.show()


