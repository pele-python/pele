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
        self.x = self.system.get_random_minimized_configuration(tol=10.).coords
    
    def do_tst(self, **kwargs):
        self.pot.nfev = 0
        opt = FindTransitionState(self.x, self.pot, **kwargs)
        ret = opt.run()
        self.assertEqual(ret.nfev, self.pot.nfev)
        self.assertTrue(ret.success)
        self.assertGreater(ret.nfev, 0)
    
    def test(self):
        self.do_tst()
 
    def test1(self):
        self.do_tst(invert_gradient=True)

    def test2(self):
        self.do_tst(lowestEigenvectorQuenchParams=dict(first_order=True)
                     )


if __name__ == "__main__":
    unittest.main()
