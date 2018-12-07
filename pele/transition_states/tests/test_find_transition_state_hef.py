import unittest
import os
import numpy as np

from pele.systems import LJCluster
from pele.transition_states import FindTransitionState
from pele.utils.xyz import read_xyz
from pele.potentials import BasePotential

class PotWrapper(BasePotential):
    def __init__(self, pot):
        self.nfev = 0
        self.pot = pot
    def getEnergy(self, coords):
        self.nfev += 1
        return self.pot.getEnergy(coords)
    def getEnergyGradient(self, coords):
        self.nfev += 1
        return self.pot.getEnergyGradient(coords)

class TestFindTransitionState(unittest.TestCase):
    def setUp(self):
        self.system = LJCluster(18)
        self.pot = PotWrapper(self.system.get_potential())
    
    def make_dimer(self, x, **kwargs):
        return FindTransitionState(x, self.pot,
                                orthogZeroEigs=self.system.get_orthogonalize_to_zero_eigenvectors(),
                                **kwargs
#                                verbosity=10
                                )
    
    def test_params(self):
        params = FindTransitionState.params()
    
    def test1(self):
        x = self.system.get_random_configuration()
        dimer = self.make_dimer(x)
        res = dimer.run()
        self.assertTrue(res.success)
        self.assertEqual(res.nfev, self.pot.nfev)
        self.assertGreater(res.nfev, 1)

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

    def test_exact_diagonalization(self):
        # get the path of the file directory
        path = os.path.dirname(os.path.abspath(__file__))
        xyz = read_xyz(open(path+"/lj18_ts.xyz", "r"))
        x = xyz.coords.flatten()
        dimer = self.make_dimer(x, hessian_diagonalization=True)
        res = dimer.run()
        self.assertTrue(res.success)
        self.assertLess(res.nsteps, 5)


class TestHEF_InvertedGradient(TestFindTransitionState):
    def make_dimer(self, x, **kwargs):
        return FindTransitionState(x, self.pot,
                                orthogZeroEigs=self.system.get_orthogonalize_to_zero_eigenvectors(),
                                invert_gradient=True,
                                )

class TestFindTransitionState_NFEV(unittest.TestCase):
    def setUp(self):
        from pele.optimize.tests.test_nfev import _PotWrapper
        np.random.seed(0)
        self.system = LJCluster(18)
        self.pot = _PotWrapper(self.system.get_potential())
        self.x = self.system.get_random_minimized_configuration(tol=10.).coords
    
    def do_check(self, **kwargs):
        self.pot.nfev = 0
        opt = FindTransitionState(self.x, self.pot, **kwargs)
        ret = opt.run()
        self.assertEqual(ret.nfev, self.pot.nfev)
        self.assertTrue(ret.success)
        self.assertGreater(ret.nfev, 0)
    
    def test(self):
        self.do_check()
 
#     def test1(self):
#         self.do_check(invert_gradient=True)

    def test2(self):
        self.do_check(lowestEigenvectorQuenchParams=dict(first_order=True)
                     )


if __name__ == "__main__":
    unittest.main()

