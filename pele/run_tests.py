import unittest

#from pele.mindist.aamindist import aaDistTest
from pele.mindist.tests import *
#from pele.mindist.minpermdist_rbmol import TestMinPermDistRBMol_OTP
#from pele.mindist.permutational_alignment import TestMinDistUtils
from pele.potentials.ATLJ import TestATLJ
from pele.potentials.lj import LJTest
from pele.potentials.ljcut import LJCutTest
from pele.landscape._graph import TestGraph
from pele.landscape._distance_graph import TestDistanceGraph
from pele.transition_states._orthogopt import TestOrthogopt
from pele.utils.hessian import TestEig
from pele.accept_tests.tests import *
from pele.storage.tests import *
from pele._test_basinhopping import TestBasinhopping

unittest.main()