import unittest
from unittest import mock
import os
import hashlib
import json
import platform
import sys
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
        # Temporary PR diagnostic: observe the original search without changing it.
        from pele.transition_states import FindLowestEigenVector
        from pele.utils import rotations

        def trace(event, **data):
            print("HEF_TRACE " + json.dumps(dict(test=self.id(), event=event, **data)), flush=True)

        def digest(array):
            return hashlib.sha256(array.tobytes()).hexdigest()

        def summary(result):
            data = {}
            for name in ("success", "nsteps", "nfev", "rms", "eigenval", "message"):
                value = getattr(result, name, None)
                if isinstance(value, np.generic):
                    value = value.item()
                if isinstance(value, float) and not np.isfinite(value):
                    value = None
                data[name] = value
            return data

        try:
            from threadpoolctl import threadpool_info
            pools = threadpool_info()
        except ImportError:
            pools = None
        state = np.random.get_state()
        rng_hash = hashlib.sha256(
            state[1].tobytes() + repr((state[0], state[2], state[3], state[4])).encode()
        ).hexdigest()
        trace("setup", python=sys.version, numpy=np.__version__, platform=platform.platform(),
              machine=platform.machine(), threadpools=pools,
              environment={name: os.environ.get(name) for name in
                           ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")},
              coords_hash=digest(self.x), coords=self.x.tolist(), rng_hash=rng_hash,
              kwargs=kwargs)

        phase = "search"
        original_run = FindLowestEigenVector.run
        original_random = rotations.vec_random_ndim

        def traced_random(shape):
            vector = original_random(shape)
            trace("random_axis", phase=phase, vector_hash=digest(vector), values=vector.tolist())
            return vector

        def traced_run(optimizer, *args, **run_kwargs):
            orthogonalizer = optimizer.eigpot.orthogZeroEigs
            label = "secondary" if getattr(orthogonalizer, "__name__", None) == "project" else "primary"
            budget = args[0] if args else run_kwargs.get("niter", optimizer.minimizer.nsteps)
            try:
                result = original_run(optimizer, *args, **run_kwargs)
            except Exception as error:
                trace("rr_error", phase=phase, label=label, budget=budget,
                      error=type(error).__name__, detail=str(error))
                raise
            trace("rr", phase=phase, label=label, budget=budget,
                  first_order=optimizer.eigpot.first_order, dx=optimizer.eigpot.diff,
                  **summary(result))
            return result

        self.pot.nfev = 0
        opt = FindTransitionState(self.x, self.pot, **kwargs)
        with mock.patch.object(rotations, "vec_random_ndim", side_effect=traced_random), \
                mock.patch.object(FindLowestEigenVector, "run", new=traced_run):
            ret = opt.run()
            trace("final", energy=float(ret.energy), coords=ret.coords.tolist(),
                  eigenvec=ret.eigenvec.tolist(), counted_nfev=self.pot.nfev, **summary(ret))
            try:
                raw_potential = self.pot.pot
                xyz = ret.coords.reshape(-1, 3)
                xyz = xyz - xyz.mean(axis=0)
                rigid = np.column_stack(
                    [np.tile(axis, len(xyz)) for axis in np.eye(3)]
                    + [np.cross(axis, xyz).ravel() for axis in np.eye(3)]
                )
                rank = int(np.linalg.matrix_rank(rigid))
                basis = np.linalg.qr(rigid, mode="complete")[0][:, rank:]
                hessian = raw_potential.getHessian(ret.coords)
                values = np.linalg.eigvalsh(basis.T @ hessian @ basis)
                trace("physical_spectrum", rigid_rank=rank, values=values.tolist(),
                      negative_modes=int(np.count_nonzero(values < -1e-5)))
                if not ret.success:
                    for budget in (100, 500):
                        phase = "fixed_replay_%d" % budget
                        params = opt.lowestEigenvectorQuenchParams.copy()
                        params["nsteps"] = budget
                        replay = FindTransitionState(
                            ret.coords.copy(), raw_potential, eigenvec0=ret.eigenvec.copy(),
                            orthogZeroEigs=opt.orthogZeroEigs,
                            lowestEigenvectorQuenchParams=params, verbosity=0,
                        )
                        result = replay._get_lowest_transverse_eigenvector(ret.coords.copy(), ret.grad.copy())
                        trace("fixed_replay", budget=budget, **summary(result))
            except Exception as error:
                trace("posthoc_error", error=type(error).__name__, detail=str(error))
        self.assertEqual(ret.nfev, self.pot.nfev)
        self.assertTrue(ret.success)
        self.assertGreater(ret.nfev, 0)

    def test(self):
        self.do_check()

    #     def test1(self):
    #         self.do_check(invert_gradient=True)

    def test2(self):
        self.do_check(lowestEigenvectorQuenchParams=dict(first_order=True))


if __name__ == "__main__":
    unittest.main()
