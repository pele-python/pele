"""Temporary observer for the installed Intel wheel; no solver settings change."""
import json
import platform
import sys
from pathlib import Path
from unittest import mock

import numpy as np
import pytest
import scipy
from threadpoolctl import threadpool_info

import pele
from pele.transition_states.tests import test_find_transition_state_hef as tests

output = open(sys.argv[1], 'w')
original_search = tests.FindTransitionState
active_case = None


def record(kind, **values):
    output.write(json.dumps(dict(kind=kind, case=active_case, **values),
                            default=lambda value: value.item() if isinstance(value, np.generic) else str(value)) + '\n')
    output.flush()


class TraceSearch(original_search):
    def _getLowestEigenVector(self, coords, iteration, gradient=None):
        overlap = super()._getLowestEigenVector(coords, iteration, gradient=gradient)
        energy, raw_gradient, hessian = self.pot.pot.getEnergyGradientHessian(coords)
        centered = coords.reshape(-1, 3) - coords.reshape(-1, 3).mean(axis=0)
        rigid = np.column_stack([
            np.tile(axis, len(centered)) for axis in np.eye(3)
        ] + [np.cross(np.tile(axis, (len(centered), 1)), centered).ravel()
             for axis in np.eye(3)])
        basis = np.linalg.qr(rigid, mode='complete')[0][:, np.linalg.matrix_rank(rigid):]
        spectrum = np.linalg.eigvalsh(basis.T @ hessian @ basis)
        curvature = float(self.eigenvec @ hessian @ self.eigenvec)
        record('primary', iteration=iteration, npositive=self.npositive,
               coords=coords.tolist(), eigenvec=self.eigenvec.tolist(),
               energy=energy, raw_rms=float(np.linalg.norm(raw_gradient) / np.sqrt(coords.size)),
               primary={key: self.leig_result.get(key) for key in ('success', 'nsteps', 'rms', 'eigenval')},
               primary_curvature=curvature,
               primary_residual=float(np.linalg.norm(basis.T @ (hessian @ self.eigenvec - curvature * self.eigenvec))),
               physical_spectrum=spectrum.tolist())
        return overlap

    def run(self):
        record('start', coords=self.coords.tolist())
        result = super().run()
        record('result', success=result.success, energy=result.energy, rms=result.rms,
               nsteps=result.nsteps, nfev=result.nfev, counted_nfev=self.pot.nfev,
               coords=result.coords.tolist(), eigenval=result.eigenval, message=result.message)
        assert result.nfev == self.pot.nfev
        return result


class TraceAccountingTests:
    def pytest_runtest_setup(self, item):
        global active_case
        if item.cls is tests.TestFindTransitionState_NFEV:
            active_case = item.nodeid
            tests.FindTransitionState = TraceSearch

    def pytest_runtest_teardown(self, item):
        tests.FindTransitionState = original_search


record('runtime', python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
       platform=platform.platform(), threads=threadpool_info())
package = Path(pele.__file__).parent
suite_status = pytest.main([str(package), '--ignore=' + str(package / 'gui'), '-q'],
                           plugins=[TraceAccountingTests()])
record('suite', status=int(suite_status))
for trial in range(10):
    for variant in ('test', 'test2'):
        active_case = f'repeat-{trial}-{variant}'
        case = tests.TestFindTransitionState_NFEV(variant)
        case.setUp()
        with mock.patch.object(tests, 'FindTransitionState', TraceSearch):
            try:
                getattr(case, variant)()
            except AssertionError as error:
                record('assertion', message=str(error))
output.close()
