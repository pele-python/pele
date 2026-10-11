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
        return FindTransitionState(
            x,
            self.pot,
            orthogZeroEigs=self.system.get_orthogonalize_to_zero_eigenvectors(),
            **kwargs
            #                                verbosity=10
        )

    def test_params(self):
        params = FindTransitionState.params()

    def test1(self):
        path = os.path.join(os.path.dirname(__file__), "lj18_ts.xyz")
        with open(path) as source:
            x = read_xyz(source).coords.ravel()
        perturbation = np.sin(np.arange(x.size) + 1.0)
        perturbation = self.system.get_orthogonalize_to_zero_eigenvectors()(perturbation, x)
        x = x + 0.01 * perturbation / np.linalg.norm(perturbation)
        dimer = self.make_dimer(x)
        res = dimer.run()
        self.assertTrue(res.success)
        self.assertEqual(res.nfev, self.pot.nfev)
        self.assertGreater(res.nfev, 1)

        # Check physical curvature independently of the iterative index validation.
        _, _, hessian = self.pot.pot.getEnergyGradientHessian(res.coords)
        centered = res.coords.reshape(-1, 3)
        centered = centered - centered.mean(axis=0)
        rigid = np.column_stack([
            np.tile(axis, len(centered)) for axis in np.eye(3)
        ] + [np.cross(np.tile(axis, (len(centered), 1)), centered).ravel()
             for axis in np.eye(3)])
        self.assertEqual(np.linalg.matrix_rank(rigid), 6)
        basis = np.linalg.qr(rigid, mode="complete")[0][:, 6:]
        eigenvalues = np.linalg.eigvalsh(basis.T @ hessian @ basis)
        self.assertLess(eigenvalues[0], -1e-5)
        self.assertGreater(eigenvalues[1], 1e-5)

    def test2(self):
        # get the path of the file directory
        path = os.path.dirname(os.path.abspath(__file__))
        xyz = read_xyz(open(path + "/lj18_ts.xyz"))
        x = xyz.coords.flatten()
        dimer = self.make_dimer(x)
        res = dimer.run()
        self.assertTrue(res.success)
        self.assertLess(res.nsteps, 5)

    #        print res

    def test_exact_diagonalization(self):
        # get the path of the file directory
        path = os.path.dirname(os.path.abspath(__file__))
        xyz = read_xyz(open(path + "/lj18_ts.xyz"))
        x = xyz.coords.flatten()
        dimer = self.make_dimer(x, hessian_diagonalization=True)
        res = dimer.run()
        self.assertTrue(res.success)
        self.assertLess(res.nsteps, 5)


class TestHEF_InvertedGradient(TestFindTransitionState):
    def make_dimer(self, x, **kwargs):
        return FindTransitionState(
            x,
            self.pot,
            orthogZeroEigs=self.system.get_orthogonalize_to_zero_eigenvectors(),
            invert_gradient=True,
            **kwargs,
        )


class TestFindTransitionState_NFEV(unittest.TestCase):
    def setUp(self):
        from pele.optimize.tests.test_nfev import _PotWrapper

        np.random.seed(0)
        self.system = LJCluster(18)
        self.pot = _PotWrapper(self.system.get_potential())
        self.x = self.system.get_random_minimized_configuration(tol=10.0).coords

    def do_check(self, **kwargs):
        self.pot.nfev = 0
        opt = FindTransitionState(self.x, self.pot, **kwargs)
        ret = opt.run()
        self.assertEqual(ret.nfev, self.pot.nfev)
        if not ret.success:
            import json
            import platform
            import sys

            _, raw_gradient, hessian = self.pot.pot.getEnergyGradientHessian(ret.coords)
            centered = ret.coords.reshape(-1, 3) - ret.coords.reshape(-1, 3).mean(axis=0)
            rigid = np.column_stack([
                np.tile(axis, len(centered)) for axis in np.eye(3)
            ] + [np.cross(np.tile(axis, (len(centered), 1)), centered).ravel()
                 for axis in np.eye(3)])
            rank = np.linalg.matrix_rank(rigid)
            basis = np.linalg.qr(rigid, mode="complete")[0][:, rank:]
            spectrum = np.linalg.eigvalsh(basis.T @ hessian @ basis)
            primary_curvature = float(ret.eigenvec @ hessian @ ret.eigenvec)
            print("HEF_TRACE " + json.dumps({
                "test": self.id(), "python": sys.version, "numpy": np.__version__,
                "platform": platform.platform(), "coords": ret.coords.tolist(),
                "energy": float(ret.energy), "rms": float(ret.rms),
                "raw_rms": float(np.linalg.norm(raw_gradient) / np.sqrt(ret.coords.size)),
                "nfev": int(ret.nfev), "counted_nfev": int(self.pot.nfev),
                "nsteps": int(ret.nsteps), "message": ret.message,
                "eigenvec": ret.eigenvec.tolist(), "eigenval": float(ret.eigenval),
                "primary": {key: opt.leig_result.get(key)
                            for key in ("success", "nsteps", "rms", "eigenval")},
                "primary_curvature": primary_curvature,
                "primary_residual": float(np.linalg.norm(
                    basis.T @ (hessian @ ret.eigenvec - primary_curvature * ret.eigenvec))),
                "physical_spectrum": spectrum.tolist(),
            }, default=lambda value: value.item() if isinstance(value, np.generic) else str(value)))
        self.assertTrue(ret.success, msg=repr(ret))
        self.assertGreater(ret.nfev, 0)

    def test(self):
        self.do_check()

    #     def test1(self):
    #         self.do_check(invert_gradient=True)

    def test2(self):
        self.do_check(lowestEigenvectorQuenchParams=dict(first_order=True))


if __name__ == "__main__":
    unittest.main()
