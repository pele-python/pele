import unittest
import numpy as np
import os

from pele.potentials._pythonpotential import CppPotentialWrapper
from pele.potentials import BasePotential, _lj_cpp
from pele.optimize import MODIFIED_FIRE_CPP

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

class TestMODIFIED_FIRE_CPP_PP(unittest.TestCase):
    def test_raises(self):
        with self.assertRaises(NotImplementedError):
            modified_fire = MODIFIED_FIRE_CPP(_xrand, _Raise())
            modified_fire.run()

class TestMODIFIED_FIRE_CPP(unittest.TestCase):
    
    def do_test(self, pot, **kwargs):
        modified_fire = MODIFIED_FIRE_CPP(_xrand, pot, **kwargs)
        res = modified_fire.run()
        self.assertAlmostEqual(res.energy, _emin, 4)
        self.assertTrue(res.success)
        self.assertLess(np.max(np.abs(res.coords - _xmin)), 1e-2)
        self.assertGreater(res.nfev, 0)
        
    def test_E(self):
        self.do_test(_E())

    def test_EG(self):
        self.do_test(_EG())

    def assert_same(self, res1, res2):
        self.assertEqual(res1.energy, res2.energy)
        self.assertEqual(res1.rms, res2.rms)
        self.assertEqual(res1.nfev, res2.nfev)
        
    def test_run_niter(self):
        modified_fire1 = MODIFIED_FIRE_CPP(_xrand, _EG())
        res1 = modified_fire1.run()
        modified_fire2 = MODIFIED_FIRE_CPP(_xrand, _EG())
        res2 = modified_fire2.run(res1.nsteps)
        self.assert_same(res1, res2)
        
    def test_run_niter2(self):
        modified_fire1 = MODIFIED_FIRE_CPP(_xrand, _EG())
        res1 = modified_fire1.run()
        modified_fire2 = MODIFIED_FIRE_CPP(_xrand, _EG())
        res2 = modified_fire2.run(res1.nsteps / 2)
        res2 = modified_fire2.run()
        self.assert_same(res1, res2)
        
    def test_run_niter3(self):
        modified_fire1 = MODIFIED_FIRE_CPP(_xrand, _EG())
        res1 = modified_fire1.run(10)
        modified_fire2 = MODIFIED_FIRE_CPP(_xrand, _EG())
        res2 = modified_fire2.run(5)
        res2 = modified_fire2.run(5)
        self.assert_same(res1, res2)
        
    def test_event_raise(self):
        class EventException(BaseException): pass
        def myevent(*args, **kwargs): raise EventException
        with self.assertRaises(EventException):
            modified_fire = MODIFIED_FIRE_CPP(_xrand, _EG(), events=[myevent])
            modified_fire.run()
    
    def test_event(self):
        self.event_called = False
        def myevent(*args, **kwargs): 
            self.event_called = True
        self.do_test(_EG(), events=[myevent])
        self.assertTrue(self.event_called)
        
class TestMODIFIED_FIRE_CPP_Raises(unittest.TestCase):
    def test_raises(self):
        pot = _lj_cpp._ErrorPotential()
        with self.assertRaises(RuntimeError):
            modified_fire = MODIFIED_FIRE_CPP(_xrand, pot)
            modified_fire.run()

#class TestMODIFIED_FIRE_CPP_PassGrad(unittest.TestCase):
#    def do_test(self, pot):
#        e, grad = pot.getEnergyGradient(_xrand)
#        modified_fire = MODIFIED_FIRE_CPP(_xrand, pot, energy=e, gradient=grad)
#        res = modified_fire.run()
#        self.assertAlmostEqual(res.energy, _emin, 4)
#        self.assertTrue(res.success)
#        self.assertLess(np.max(np.abs(res.coords - _xmin)), 1e-2)
#        self.assertGreater(res.nfev, 0)
#        
#    def test_E(self):
#        self.do_test(_E())
#
#    def test_EG(self):
#        self.do_test(_EG())

if __name__ == "__main__":
    unittest.main()
