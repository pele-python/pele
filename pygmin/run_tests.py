import unittest

from pygmin.mindist.aamindist import aaDistTest
from pygmin.mindist.minpermdist_stochastic import TestMinPermDistStochastic_BLJ
from pygmin.mindist.minpermdist_rbmol import TestMinPermDistRBMol_OTP
from pygmin.mindist.permutational_alignment import TestMinDistUtils
from pygmin.potentials.ATLJ import TestATLJ
from pygmin.potentials.lj import LJTest
from pygmin.potentials.ljcut import LJCutTest
from pygmin.landscape._graph import TestGraph
from pygmin.landscape._distance_graph import TestDistanceGraph
from pygmin.transition_states._orthogopt import TestOrthogopt

unittest.main()