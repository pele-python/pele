import numpy as np
import unittest
import logging



class _BaseTest(unittest.TestCase):
    # def energy_test(self, x, e):
    # e = self.pot.getEnergy(x)
    # print e
    #        self.assertAlmostEqual(e, self.target_E, 4)

    def grad_t(self, x):
        log = logging.getLogger("BaseTest.grad_t")
        e, g = self.pot.getEnergyGradient(x)
        e1 = self.pot.getEnergy(x)
        numerical_g = self.pot.NumericalDerivative(x)
        log.debug("g= %r", g)
        log.debug("numerical_g= %r", numerical_g)
        self.assertLess(np.max(np.abs(g - numerical_g)), 1e-3)
        self.assertAlmostEqual(e, e1, 4)

    def test_e_min(self):
        log = logging.getLogger("BaseTest.test_e_min")
        e = self.pot.getEnergy(self.xmin)
        log.debug("e= %r", e)
        self.assertAlmostEqual(e, self.Emin, 4)

    def test_grad_min(self):
        log = logging.getLogger("BaseTest.test_gra_min")
        e, g = self.pot.getEnergyGradient(self.xmin)
        log.debug("e= %r", e)
        log.debug("g= %r", g)
        self.assertAlmostEqual(e, self.Emin, 4)
        self.assertLess(np.max(np.abs(g)), 1e-3)
        self.grad_t(self.xmin)

    def test_hess_min(self):
        log = logging.getLogger("BaseTest.test_hess_min")
        h = self.pot.getHessian(self.xmin)
        eigenvals = np.linalg.eigvals(h)
        log.debug("e= %r", eigenvals)
        self.assertGreater(np.min(eigenvals), -1e-4)

    def test_hess_analytical_against_numerical(self):
        log = logging.getLogger("BaseTest.test_hess_analytical_against_numerical")
        h = self.pot.getHessian(self.xmin)
        h_num = self.pot.NumericalHessian(self.xmin)
        h = h.reshape(-1).copy()
        h_num = h_num.reshape(-1).copy()
        np.allclose(h, h_num, rtol=1e-8)

    def test_random(self):
        self.grad_t(self.xmin + self.xrandom)


def assert_arrays_almost_equal(self, v1, v2, **kwargs):
    self.assertEqual(v1.shape, v2.shape)
    for x, y in zip(v1.reshape(-1), v2.reshape(-1)):
        self.assertAlmostEqual(x, y, **kwargs)


class _TestConfiguration(unittest.TestCase):
    """test a potential at a single configuration single configuration
    
    In setUp() you should define attributes self.pot, self.x0 and the energy of that position x0
    """
    pot = None
    x0 = None
    e0 = None
    ae_kwargs = dict(places=3)  # kwargs passed to assertAlmostEqual

    def compare_arrays(self, v1, v2):
        assert_arrays_almost_equal(self, v1, v2, **self.ae_kwargs)

    def test_energy(self):
        self.assertAlmostEqual(self.pot.getEnergy(self.x0), self.e0, **self.ae_kwargs)

    def test_energy_gradient(self):
        e, g = self.pot.getEnergyGradient(self.x0)
        self.assertAlmostEqual(e, self.e0, **self.ae_kwargs)
        gnum = self.pot.NumericalDerivative(self.x0)
        self.compare_arrays(g, gnum)

    def test_energy_gradient_hessian(self):
        e, g, hess = self.pot.getEnergyGradientHessian(self.x0)
        self.assertAlmostEqual(e, self.e0, **self.ae_kwargs)
        gnum = self.pot.NumericalDerivative(self.x0)
        self.compare_arrays(g, gnum)
        hnum = self.pot.NumericalHessian(self.x0)
        self.compare_arrays(hess.reshape(-1), hnum.reshape(-1))

