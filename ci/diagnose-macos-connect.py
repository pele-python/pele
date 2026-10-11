"""Trace the seeded connection test without changing the installed wheel."""
import hashlib
import json
import logging
import platform
import sys
import traceback

import numpy as np
import scipy

import pele.landscape.local_connect as local_connect
from pele.mindist.minpermdist_stochastic import MinPermDistCluster
from pele.systems import LJCluster
from pele.optimize import LBFGS
from pele.transition_states._transverse_walker import _TransversePotential
from pele.transition_states._transition_state_refinement import FindTransitionState
from pele.transition_states import NEBDriver
from pele.landscape.tests.test_double_ended_connect import TestDoubleEndedConnect
from pele.transition_states.tests.test_NEB import _x1, _x2

logging.basicConfig(level=logging.INFO)
trace = open("connect-trace.jsonl", "w", encoding="utf-8")


def record(kind, **data):
    trace.write(json.dumps({"kind": kind, **data}) + "\n")
    trace.flush()
    print(kind, {key: value for key, value in data.items()
                 if key not in ("coords", "path", "energies", "eigenvector", "displacement")}, flush=True)


record("runtime", python=sys.version, machine=platform.machine(),
       numpy=np.__version__, scipy=scipy.__version__)
np.show_config()
np.show_runtime()
scipy.show_config()
pot = LJCluster(13).get_potential()
for label, coords in (("x1", _x1), ("x2", _x2)):
    energy, gradient, hessian = pot.getEnergyGradientHessian(coords)
    record("endpoint", label=label, energy=float(energy),
           rms=float(np.linalg.norm(gradient) / np.sqrt(coords.size)),
           eigenvalues=np.linalg.eigvalsh(hessian).tolist())

align = MinPermDistCluster.__call__


def traced_align(self, *args, **kwargs):
    result = align(self, *args, **kwargs)
    record("alignment", distance=float(result[0]), coords=[x.tolist() for x in result[1:]])
    return result


MinPermDistCluster.__call__ = traced_align
run_neb = NEBDriver.run
band = 0


def traced_neb(self):
    global band
    band += 1
    number = band

    def snapshot(path=None, energies=None, stepnum=None, rms=None, event=None, **kwargs):
        if event != "update" or stepnum in (1, 2, 10, 100, 500, 1000):
            record("band", number=number, event=event, evaluation=stepnum,
                   rms=float(rms), digest=hashlib.sha256(path.tobytes()).hexdigest(),
                   path=path.tolist(), energies=energies.tolist())

    self.update_event.connect(snapshot)
    record("neb_parameters", number=number, quench=self.quenchParams,
           iter_density=self.iter_density, factor=self.factor)
    return run_neb(self)


NEBDriver.run = traced_neb
refine = local_connect._refineTS
candidate = 0
tracing_candidate = False
adjust_step = LBFGS.adjustStepSize
get_eigenvector = FindTransitionState._getLowestEigenVector
transverse_steps = []
eigenvectors = []
saddles = []


def traced_adjust_step(self, coords, energy, gradient, step):
    before = coords.copy()
    result = adjust_step(self, coords, energy, gradient, step)
    if tracing_candidate and isinstance(self.pot, _TransversePotential):
        vector = self.pot.vector
        displacement = result[0] - before
        length = float(np.linalg.norm(displacement))
        transverse_steps.append((before, vector.copy(), displacement.copy(), gradient.copy()))
        record("transverse_displacement", iteration=self.iter_number,
               length=length, parallel=float(np.dot(displacement, vector)),
               relative_parallel=float(np.dot(displacement, vector)) / max(length, 1e-300),
               gradient_parallel=float(np.dot(gradient, vector)),
               energy_change=float(result[1] - energy), coords=before.tolist(),
               eigenvector=vector.tolist(), displacement=displacement.tolist())
    return result


def traced_eigenvector(self, coords, iteration, gradient=None):
    result = get_eigenvector(self, coords, iteration, gradient=gradient)
    if tracing_candidate:
        eigenvectors.append((iteration, coords.copy(), self.eigenval, self.eigenvec.copy()))
        record("lowest_eigenvector", iteration=iteration, coords=coords.tolist(),
               eigenvalue=float(self.eigenval), eigenvector=self.eigenvec.tolist(),
               success=bool(self.leig_result.success), rms=float(self.leig_result.rms))
    return result


LBFGS.adjustStepSize = traced_adjust_step
FindTransitionState._getLowestEigenVector = traced_eigenvector


def traced_refine(*args, **kwargs):
    global candidate, tracing_candidate
    candidate += 1
    tracing_candidate = candidate == 1
    try:
        result = refine(*args, **kwargs)
    finally:
        tracing_candidate = False
    for label, ret in zip(("ts", "minimum1", "minimum2"), result[1:]):
        if ret is not None:
            if label == "ts":
                saddles.append(ret.coords.copy())
            record("refinement", label=label, success=bool(ret.success),
                   energy=float(ret.energy), rms=float(ret.rms), coords=ret.coords.tolist())
    return result


local_connect._refineTS = traced_refine
try:
    TestDoubleEndedConnect("test2").test2()
except Exception:
    record("result", success=False, error=traceback.format_exc())
    raise
else:
    record("result", success=True)
finally:
    # Evaluate physical checks afterward so they cannot change the search trajectory.
    for step, (coords, vector, displacement, projected_gradient) in enumerate(transverse_steps):
        _, true_gradient = pot.getEnergyGradient(coords)
        record("physical_slope", step=step,
               true_slope=float(np.dot(true_gradient, displacement)),
               projected_slope=float(np.dot(projected_gradient, displacement)),
               excluded_slope=float(np.dot(true_gradient, vector) * np.dot(vector, displacement)))
    for label, points in (("eigenpair", eigenvectors), ("saddle", enumerate(saddles))):
        for item in points:
            iteration, coords = item[:2]
            _, gradient, hessian = pot.getEnergyGradientHessian(coords)
            centered = coords.reshape(-1, 3) - coords.reshape(-1, 3).mean(axis=0)
            rigid = np.column_stack([
                np.tile(axis, len(centered)) for axis in np.eye(3)
            ] + [np.cross(np.tile(axis, (len(centered), 1)), centered).ravel()
                 for axis in np.eye(3)])
            assert np.linalg.matrix_rank(rigid) == 6
            basis = np.linalg.qr(rigid, mode="complete")[0][:, 6:]
            assert np.allclose(basis.T @ rigid, 0, atol=1e-12)
            projected_hessian = basis.T @ hessian @ basis
            values = np.linalg.eigvalsh(projected_hessian)
            details = dict(iteration=iteration, eigenvalues=values.tolist(),
                           negative_modes=int(np.count_nonzero(values < -1e-5)),
                           physical_gradient_rms=float(np.linalg.norm(gradient) / np.sqrt(coords.size)))
            if label == "eigenpair":
                _, _, estimate, vector = item
                vector = basis.T @ vector
                vector /= np.linalg.norm(vector)
                rayleigh = float(vector @ projected_hessian @ vector)
                details.update(estimate=float(estimate), rayleigh=rayleigh,
                               residual=float(np.linalg.norm(projected_hessian @ vector - rayleigh * vector)),
                               lowest_error=float(rayleigh - values[0]))
            record("exact_" + label, **details)
    trace.close()
