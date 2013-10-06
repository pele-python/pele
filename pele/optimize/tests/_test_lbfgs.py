import unittest

from pele.optimize import LBFGS, MYLBFGS
from pele.systems import LJCluster

class TestLBFGS_State(unittest.TestCase):
    def setUp(self):
        self.system = LJCluster(13)
        self.x = self.system.get_random_configuration()
        self.pot = self.system.get_potential()
        self.minimizer = LBFGS(self.x, self.pot)
    
    def test_state(self):
        # do several minimization iterations
        for i in xrange(10):
            self.minimizer.one_iteration()
        
        # get the state and save it
        ret = self.minimizer.get_result()
        state = self.minimizer.get_state()
        x1 = ret.coords.copy()
        
        # do several more iteration steps
        for i in xrange(10):
            self.minimizer.one_iteration()
        
        # now make a new minimizer and do several iterations
        minimizer2 = LBFGS(x1, self.pot)
        minimizer2.set_state(state)
        for i in xrange(10):
            minimizer2.one_iteration()
        
        # test that the two minimizers are in the same state
        ret1 = self.minimizer.get_result()
        ret2 = minimizer2.get_result()
        self.assertEqual(ret1.energy, ret2.energy)
        self.assertTrue((ret1.coords == ret2.coords).all())
        
        state1 = self.minimizer.get_state()
        state2 = minimizer2.get_state()
        
        self.assertTrue((state1.y == state2.y).all())
        self.assertTrue((state1.s == state2.s).all())
        self.assertTrue((state1.rho == state2.rho).all())
        self.assertTrue((state1.Xold == state2.Xold).all())
        self.assertTrue((state1.Gold == state2.Gold).all())
        self.assertEqual(state1.H0, state2.H0)
        self.assertEqual(state1.k, state2.k)
        
#class TestLBFGS_wolfe(unittest.TestCase):
#    def setUp(self):
#        self.system = LJCluster(13)
#        self.x = self.system.get_random_configuration()
#        self.pot = self.system.get_potential()
#    
#    def test(self):
#        minimizer = LBFGS(self.x.copy(), self.pot, wolfe=True, debug=True)
#        ret = minimizer.run()
#        self.assertTrue(ret.success)
#        
#        print "\n\n"
#        minimizer = LBFGS(self.x.copy(), self.pot, wolfe=False, debug=True)
#        ret_nowolfe = minimizer.run()
#        self.assertTrue(ret_nowolfe.success)
#        
#        print "nfev wolfe, nowolfe", ret.nfev, ret_nowolfe.nfev, ret.energy, ret_nowolfe.energy
  

class TestLBFGSCython(unittest.TestCase):
    def setUp(self):
        self.system = LJCluster(13)
        self.x = self.system.get_random_configuration()
        self.pot = self.system.get_potential()
    
    def test(self):
        minimizer = LBFGS(self.x.copy(), self.pot, cython=True, debug=True)
        ret = minimizer.run()
        m2 = LBFGS(self.x.copy(), self.pot, cython=False, debug=True)
        ret2 = m2.run()
        
        print "cython", ret.nfev, ret2.nfev
        self.assertEqual(ret.nfev, ret2.nfev)
        self.assertAlmostEqual(ret.energy, ret2.energy, 5)

if __name__ == "__main__":
    unittest.main()
        
        
        