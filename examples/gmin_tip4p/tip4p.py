from math import sin, cos, pi
import numpy as np
from pygmin.potentials import GMINPotential
from pygmin.angleaxis import RigidFragment, RBSystem
import gmin_ as GMIN
from pygmin.utils import rotations

def swap(x1, x2):
    tmp = x1.copy()
    x1[:] = x2
    x2[:] = tmp

GMIN.initialize()

pot = GMINPotential(GMIN)
coords = pot.getCoords()

nrigid = coords.size / 6
print nrigid
water = RigidFragment()

# define water molecule
rho   = 0.9572
theta = 104.52/180.0*pi      
water.add_atom("O", np.array([0., 0., 0.]), 16.)
water.add_atom("H", rho*np.array([0.0, sin(0.5*theta), cos(0.5*theta)]), 1.)
water.add_atom("H", rho*np.array([0.0, -sin(0.5*theta), cos(0.5*theta)]), 1.)
water.finalize_setup()

system = RBSystem()
system.add_sites([water for i in xrange(nrigid)])

ca = system.coords_adapter(coords)

print pot.getEnergy(coords)
print ca.rotRigid[-1], coords[-3:]
#ca.rotRigid[1] = rotations.rotate_aa(ca.rotRigid[1], np.array([pi, 0., 0.]))
coords[-3:] = rotations.rotate_aa(coords[-3:], np.array([pi, 0., 0.]))
#ca.rotRigid[2] = rotations.rotate_aa(np.array([pi, 0., 0.]), ca.rotRigid[2])
#swap(ca.posRigid[0], ca.posRigid[1])
#swap(ca.rotRigid[0], ca.rotRigid[1])
print pot.getEnergy(coords)
#ca.rotRigid[2] = rotations.rotate_aa(ca.rotRigid[2], np.array([pi, 0., 0.]))
#print pot.getEnergy(coords)
