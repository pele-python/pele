import numpy as np
import copy

from . import CoordsAdapter
from pygmin.utils import rotations
from pygmin.NEB import interpolate_linear

def interpolate_angleaxis(initial, final, t):
    conf = initial.copy()
    for i in xrange(conf.shape[0]):
        conf[i] = rotations.q2aa(rotations.q_slerp(rotations.aa2q(initial[i]),
                                rotations.aa2q(final[i]), t))
    