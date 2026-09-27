import pytest
import os
import logging

import numpy as np
try:
    from unittest import mock
except ImportError:
    from unittest import mock

from pele.optimize._quench import lbfgs_cpp
from . import _base_test


def read_xyzdr(fname, bdim=3):
    coords = []
    radii = []
    rattlers = []
    f = open(fname)
    while True:
        xyzdr = f.readline()
        if not xyzdr:
            break
        x, y, z, d, r = xyzdr.split()
        coords.extend([float(x), float(y), float(z)])
        radii.extend([float(d) / 2])
        for _ in range(bdim):
            rattlers.extend([float(r)])
    return np.array(coords), np.array(radii), np.array(rattlers)


def minimize(coords, pot):
    result = lbfgs_cpp(coords, pot)
    # result = modifiedfire_cpp(coords, pot)
    return result.coords, result.energy, result.grad, result.rms


@pytest.fixture
def potential_and_minima():
    _inversepower_hs_cpp = pytest.importorskip("pele.potentials._inversepower_hs_cpp")
    current_dir = os.path.dirname(__file__)
    xyz, hs_radii, rattlers = read_xyzdr(
        current_dir + "/_hswca20_min2.xyzdr"
    )
    sigma = 0.205071132088
    boxv = np.array([6.26533756282, 6.26533756282, 6.26533756282])
    pow = 4
    eps = 1
    pot = _inversepower_hs_cpp.InversePowerHS(pow, eps, sigma, hs_radii, boxvec=boxv)
    natoms = 20
    result = minimize(xyz, pot)
    
    xmin = result[0]
    Emin = result[1]
    xrandom = np.random.uniform(-1, 1, len(xyz)) * 1e-2

    return pot, xmin, Emin, xrandom, natoms


class TestInversePowerHSCPP:
    def test_e_min(self, potential_and_minima):
        pot, xmin, Emin, _, _ = potential_and_minima
        e = pot.getEnergy(xmin)
        assert np.isclose(e, Emin, atol=1e-4)

    def test_grad_min(self, potential_and_minima):
        pot, xmin, Emin, _, _ = potential_and_minima
        e, g = pot.getEnergyGradient(xmin)
        assert np.isclose(e, Emin, atol=1e-4)
        assert np.max(np.abs(g)) < 1e-3

    def test_grad_t(self, potential_and_minima):
        pot, xmin, Emin, _, _ = potential_and_minima
        e, g = pot.getEnergyGradient(xmin)
        numerical_g = pot.NumericalDerivative(xmin)
        assert np.max(np.abs(g - numerical_g)) < 1e-3
        e1 = pot.getEnergy(xmin)
        assert np.isclose(e, e1, atol=1e-4)

    def test_hess_min(self, potential_and_minima):
        pot, xmin, _, _, _ = potential_and_minima
        e, g, h = pot.getEnergyGradientHessian(xmin)
        eigenvals = np.linalg.eigvals(h)
        assert np.min(eigenvals) > -1e-4

    def test_hess_analytical_against_numerical(self, potential_and_minima):
        pot, xmin, _, _, _ = potential_and_minima
        e, g, h = pot.getEnergyGradientHessian(xmin)
        h_num = pot.NumericalHessian(xmin)
        h = h.reshape(-1).copy()
        h_num = h_num.reshape(-1).copy()
        assert np.allclose(h, h_num, rtol=1e-8)

    def test_random(self, potential_and_minima):
        pot, xmin, _, xrandom, _ = potential_and_minima
        x = xmin + xrandom
        e, g = pot.getEnergyGradient(x)
        numerical_g = pot.NumericalDerivative(x)
        assert np.max(np.abs(g - numerical_g)) < 1e-3
        e1 = pot.getEnergy(x)
        assert np.isclose(e, e1, atol=1e-4)


if __name__ == "__main__":
    logging.basicConfig(filename="hs_wca_cpp.log", level=logging.DEBUG)
    # This part is for running with unittest, which is not the case anymore.
    # To run with pytest, just run `pytest` in the terminal.
    pass 