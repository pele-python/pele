import unittest
import os

from pele.systems import LJCluster
from pele.transition_states import FindTransitionState
from pele.utils.xyz import read_xyz

class TestFindTransitionState(unittest.TestCase):
    def setUp(self):
        self.system = LJCluster(18)
        self.pot = self.system.get_potential()
    
    def make_dimer(self, x):
        return FindTransitionState(x, self.pot,
                                orthogZeroEigs=self.system.get_orthogonalize_to_zero_eigenvectors(),
#                                verbosity=10
                                )
    
    def test_params(self):
        params = FindTransitionState.params()
    
    def test1(self):
        x = self.system.get_random_configuration()
        dimer = self.make_dimer(x)
        res = dimer.run()
        self.assertTrue(res.success)

    def test2(self):
        # get the path of the file directory
        path = os.path.dirname(os.path.abspath(__file__))
        xyz = read_xyz(open(path+"/lj18_ts.xyz", "r"))
        x = xyz.coords.flatten()
        dimer = self.make_dimer(x)
        res = dimer.run()
        self.assertTrue(res.success)
        self.assertLess(res.nsteps, 5)
#        print res

class TestHEF_InvertedGradient(TestFindTransitionState):
    def make_dimer(self, x):
        return FindTransitionState(x, self.pot,
                                orthogZeroEigs=self.system.get_orthogonalize_to_zero_eigenvectors(),
                                invert_gradient=True,
                                )

class TestFindTransitionState_NFEV(unittest.TestCase):
    def setUp(self):
        from pele.optimize.tests._test_nfev import _PotWrapper
        self.system = LJCluster(18)
        self.pot = _PotWrapper(self.system.get_potential())
    
    def test(self):
        x = self.system.get_random_configuration()
        opt = FindTransitionState(x, self.pot)
        ret = opt.run()
        self.assertEqual(ret.nfev, self.pot.nfev)
        self.assertTrue(ret.success)
        self.assertGreater(ret.nfev, 0)
 


if __name__ == "__main__":
    unittest.main()
