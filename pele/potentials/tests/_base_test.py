import numpy as np
import unittest
import logging
from itertools import izip

class _BaseTest(unittest.TestCase):
    
#    def energy_test(self, x, e):
#        e = self.pot.getEnergy(x)
#        print e
#        self.assertAlmostEqual(e, self.target_E, 4)
    
    def grad_t(self, x):
        log= logging.getLogger( "BaseTest.grad_t" )
        e, g = self.pot.getEnergyGradient(x)
        e1 = self.pot.getEnergy(x)
        numerical_g = self.pot.NumericalDerivative(x)
        log.debug( "g= %r", g )
        log.debug( "numerical_g= %r", numerical_g )
        self.assertLess(np.max(np.abs(g - numerical_g)), 1e-3)
        self.assertAlmostEqual(e, e1, 4)
    
    def test_e_min(self):
        log= logging.getLogger( "BaseTest.test_e_min" )
        e = self.pot.getEnergy(self.xmin)
        log.debug( "e= %r", e )
        self.assertAlmostEqual(e, self.Emin, 4)
        
    def test_grad_min(self):
        log= logging.getLogger( "BaseTest.test_gra_min" )
        e, g = self.pot.getEnergyGradient(self.xmin)
        log.debug( "e= %r", e )
        log.debug( "g= %r", g)
        self.assertAlmostEqual(e, self.Emin, 4)
        self.assertLess(np.max(np.abs(g)), 1e-3)
        self.grad_t(self.xmin)

    def test_hess_min(self):
        log= logging.getLogger( "BaseTest.test_hess_min" )
        h = self.pot.getHessian(self.xmin)
        eigenvals = np.linalg.eigvals(h)
        log.debug( "e= %r", eigenvals )
        self.assertGreater(np.min(eigenvals), -1e-4)
    
    def test_hess_analytical_against_numerical(self):
        log= logging.getLogger( "BaseTest.test_hess_analytical_against_numerical" )
        h = self.pot.getHessian(self.xmin)
        h_num = self.pot.NumericalHessian(self.xmin)
        h = h.reshape(-1).copy()
        h_num = h_num.reshape(-1).copy()
        np.allclose(h,h_num,rtol=1e-8)
    
    def test_random(self):
        self.grad_t(self.xmin+self.xrandom)
    

    
