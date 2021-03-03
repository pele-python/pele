from __future__ import absolute_import
from math import sin, cos, pi
import numpy as np
from .rigidbody import RigidFragment

def create_water():
    water = RigidFragment()
    rho   = 0.9572
    theta = 104.52/180.0*pi      
    water.add_atom("O", np.array([0., 0., 0.]), 16.)
    water.add_atom("H", rho*np.array([0.0, sin(0.5*theta), cos(0.5*theta)]), 1.)
    water.add_atom("H", rho*np.array([0.0, -sin(0.5*theta), cos(0.5*theta)]), 1.)
    water.finalize_setup()
    return water

