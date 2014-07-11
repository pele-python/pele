import unittest
import numpy as np
import os

from pele.potentials._pythonpotential import CppPotentialWrapper
from pele.potentials import BasePotential, _lj_cpp
from pele.optimize import LBFGS_CPP

ndof = 4
_xrand = np.random.uniform(-1,1,[ndof])
_xmin = np.zeros(ndof)
_emin = 0.

class _E(BasePotential):
    def getEnergy(self, x):
        return np.dot(x, x)

class _EG(object):
    def getEnergy(self, x):
        return np.dot(x, x)
    def getEnergyGradient(self, x):
        return self.getEnergy(x), 2. * x

class _Raise(BasePotential):
    def getEnergy(self, x):
        raise NotImplementedError
    def getEnergyGradient(self, x):
        raise NotImplementedError

class TestLBFGS_CPP_PP(unittest.TestCase):
    def test_raises(self):
        with self.assertRaises(NotImplementedError):
            lbfgs = LBFGS_CPP(_xrand, _Raise())
            lbfgs.run()
    
            

class TestLBFGS_CPP(unittest.TestCase):
    def do_check(self, pot, **kwargs):
        lbfgs = LBFGS_CPP(_xrand, pot, **kwargs)
        res = lbfgs.run()
        self.assertAlmostEqual(res.energy, _emin, 4)
        self.assertTrue(res.success)
        self.assertLess(np.max(np.abs(res.coords - _xmin)), 1e-2)
        self.assertGreater(res.nfev, 0)
        
    def test_E(self):
        self.do_check(_E())

    def test_EG(self):
        self.do_check(_EG())

    def assert_same(self, res1, res2):
        self.assertEqual(res1.energy, res2.energy)
        self.assertEqual(res1.rms, res2.rms)
        self.assertEqual(res1.nfev, res2.nfev)
        
    def test_run_niter(self):
        lbfgs1 = LBFGS_CPP(_xrand, _EG())
        res1 = lbfgs1.run()
        lbfgs2 = LBFGS_CPP(_xrand, _EG())
        res2 = lbfgs2.run(res1.nsteps)
        self.assert_same(res1, res2)
        
    def test_run_niter2(self):
        lbfgs1 = LBFGS_CPP(_xrand, _EG())
        res1 = lbfgs1.run()
        lbfgs2 = LBFGS_CPP(_xrand, _EG())
        res2 = lbfgs2.run(res1.nsteps / 2)
        res2 = lbfgs2.run()
        self.assert_same(res1, res2)
        
    def test_run_niter3(self):
        lbfgs1 = LBFGS_CPP(_xrand, _EG())
        res1 = lbfgs1.run(10)
        lbfgs2 = LBFGS_CPP(_xrand, _EG())
        res2 = lbfgs2.run(5)
        res2 = lbfgs2.run(5)
        self.assert_same(res1, res2)
    
    def test_result(self):
        lbfgs = LBFGS_CPP(_xrand, _EG())
        res = lbfgs.one_iteration()
        self.assertIn("H0", res)
        self.assertIn("energy", res)
        self.assertIn("grad", res)
        self.assertIn("success", res)
        self.assertIn("coords", res)
        self.assertIn("rms", res)
        self.assertIn("nsteps", res)
        self.assertIn("nfev", res)

      
    def test_event_raise(self):
        class EventException(BaseException): pass
        def myevent(*args, **kwargs): raise EventException
        with self.assertRaises(EventException):
            lbfgs = LBFGS_CPP(_xrand, _EG(), events=[myevent])
            lbfgs.run()
    
    def test_event(self):
        self.event_called = False
        def myevent(*args, **kwargs): 
            self.event_called = True
        self.do_check(_EG(), events=[myevent])
        self.assertTrue(self.event_called)
    
    def test_logger(self):
        self.do_check(_EG(), logger=True)

    def test_rel_energy(self):
        self.do_check(_EG(), rel_energy=True)


class TestLBFGS_CPP_PassGrad(unittest.TestCase):
    def do_check(self, pot):
        e, grad = pot.getEnergyGradient(_xrand)
        lbfgs = LBFGS_CPP(_xrand, pot, energy=e, gradient=grad)
        res = lbfgs.run()
        self.assertAlmostEqual(res.energy, _emin, 4)
        self.assertTrue(res.success)
        self.assertLess(np.max(np.abs(res.coords - _xmin)), 1e-2)
        self.assertGreater(res.nfev, 0)
        
    def test_E(self):
        self.do_check(_E())

    def test_EG(self):
        self.do_check(_EG())

class TestLBFGS_CPP_Raises(unittest.TestCase):
    def test_raises(self):
        pot = _lj_cpp._ErrorPotential()
        with self.assertRaises(RuntimeError):
            lbfgs = LBFGS_CPP(_xrand, pot)
            lbfgs.run()


class TestLBFGS_CPP_LJ(unittest.TestCase):
    def setUp(self):
        from pele.potentials import LJ
        self.x0 = np.zeros(6)
        self.x0[0] = 2.
        self.pot = LJ()
    
    def test_reset(self):
        lbfgs1 = LBFGS_CPP(self.x0, self.pot)
        lbfgs1.run()
        res1 = lbfgs1.get_result()
        
        x2 = self.x0.copy()
        x2[1] = 2.
        lbfgs2 = LBFGS_CPP(x2, self.pot)
        H0 = lbfgs2.get_result()["H0"]
        lbfgs2.run()
        lbfgs2.reset(self.x0)
        lbfgs2.set_H0(H0)
        lbfgs2.run()
        res2 = lbfgs2.get_result()
        
        self.assertEqual(res1.rms, res2.rms)
        self.assertEqual(res1.H0, res2.H0)
        self.assertEqual(res1.nfev, res2.nfev)
        self.assertEqual(res1.nsteps, res2.nsteps)
        self.assertTrue(np.all(res1.coords == res2.coords))

if __name__ == "__main__":
    unittest.main()
