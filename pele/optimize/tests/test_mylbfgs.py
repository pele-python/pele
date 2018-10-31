import unittest

from pele.optimize import MYLBFGS
from pele.systems import LJCluster


class TestMYLBFGS_State(unittest.TestCase):
    def setUp(self):
        self.system = LJCluster(13)
        self.x = self.system.get_random_configuration()
        self.pot = self.system.get_potential()
        self.minimizer = MYLBFGS(self.x, self.pot)

    def test_state(self):
        # do several minimization iterations
        for i in range(10):
            self.minimizer.one_iteration()

        # get the state and save it
        ret = self.minimizer.get_result()
        state = self.minimizer.get_state()
        x1 = ret.coords.copy()

        # do several more iteration steps
        for i in range(10):
            self.minimizer.one_iteration()

        # now make a new minimizer and do several iterations
        minimizer2 = MYLBFGS(x1, self.pot)
        minimizer2.set_state(state)
        for i in range(10):
            minimizer2.one_iteration()

        # test that the two minimizers are in the same state
        ret1 = self.minimizer.get_result()
        ret2 = minimizer2.get_result()
        self.assertEqual(ret1.energy, ret2.energy)
        self.assertTrue((ret1.coords == ret2.coords).all())

        state1 = self.minimizer.get_state()
        state2 = minimizer2.get_state()

        self.assertTrue((state1.W == state2.W).all())
        self.assertTrue((state1.dXold == state2.dXold).all())
        self.assertTrue((state1.dGold == state2.dGold).all())
        self.assertEqual(state1.H0, state2.H0)
        self.assertEqual(state1.point, state2.point)
        self.assertEqual(state1.iter, state2.iter)


if __name__ == "__main__":
    unittest.main()
        
        
        
