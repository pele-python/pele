import contextlib
import hashlib
import json
import os
import inspect
import platform
import textwrap
import scipy
import struct
import sys

import numpy as np
from pele.systems import LJCluster
from pele.optimize import LBFGS_CPP
from pele.optimize.tests.test_nfev import _PotWrapper
from pele.transition_states import FindTransitionState, FindLowestEigenVector
import pele.utils.rotations as rotations

out = open(sys.argv[1], 'w')


def emit(kind, **data):
    out.write(json.dumps(dict(kind=kind, **data), allow_nan=True) + '\n')
    out.flush()


def digest(array):
    return hashlib.sha256(np.asarray(array, dtype='<f8').tobytes()).hexdigest()


def state():
    name, values, position, has_gauss, cached = np.random.get_state()
    encoded = name.encode() + values.astype('<u4', copy=False).tobytes()
    encoded += struct.pack('<qqd', position, has_gauss, cached)
    return dict(hash=hashlib.sha256(encoded).hexdigest(), position=int(position),
                has_gauss=int(has_gauss), cached_gaussian=float(cached))


def result_fields(result):
    keys = ['success', 'energy', 'eigenval', 'rms', 'nsteps', 'nfev', 'message', 'H0']
    data = {}
    for key in keys:
        if key not in result:
            continue
        value = result[key]
        if isinstance(value, np.generic):
            value = value.item()
        data[key] = value
    for key in ['coords', 'eigenvec']:
        if key in result:
            data[key] = result[key].tolist()
            data[key + '_hash'] = digest(result[key])
    return data


emit('runtime', python=sys.version, numpy=np.__version__, scipy=scipy.__version__, machine=platform.machine(), platform=platform.platform(), compiler=os.environ.get('CC'))
np.show_config()
np.show_runtime()
scipy.show_config()
np.random.seed(0)
emit('seeded', rng=state())
system = LJCluster(18)
emit('system_created', rng=state())
pot = system.get_potential()
random_configuration = system.get_random_configuration


def configuration():
    replay = np.random.RandomState()
    replay.set_state(np.random.get_state())
    uniform = replay.uniform(-1, 1, 54)
    raw = random_configuration()
    emit('raw_configuration', uniform=uniform.tolist(), uniform_hash=digest(uniform),
         scale_pow=float(18) ** (1.0 / 3), scale_product=0.7 * float(18) ** (1.0 / 3),
         coords=raw.tolist(), coords_hash=digest(raw), rng=state())
    return raw


system.get_random_configuration = configuration
with open(os.devnull, 'w') as sink, contextlib.redirect_stdout(sink):
    loose = system.get_random_minimized_configuration(tol=10.0)
emit('loose_quench_default', **result_fields(loose), rng=state())
after_quench = np.random.get_state()

# Trace a separate quench, then prove that Python event delivery preserved the result.
replay = np.random.RandomState(0)
raw = replay.uniform(-1, 1, 54) * 0.7 * float(18) ** (1.0 / 3)
quench_steps = []


def quench_event(**values):
    coords = values['coords']
    quench_steps.append(coords.copy())
    emit('loose_quench_step', step=len(quench_steps), energy=float(values['energy']),
         rms=float(values['rms']), coords=coords.tolist(), coords_hash=digest(coords))


with open(os.devnull, 'w') as sink, contextlib.redirect_stdout(sink):
    observed_loose = LBFGS_CPP(raw, pot, tol=10.0, events=[quench_event]).run()
emit('loose_quench_observer_comparison', equal_coords=bool(np.array_equal(loose.coords, observed_loose.coords)),
     max_abs_difference=float(np.max(np.abs(loose.coords - observed_loose.coords))),
     equal_energy=bool(loose.energy == observed_loose.energy),
     equal_rms=bool(loose.rms == observed_loose.rms),
     equal_rng=bool(state()['hash'] == hashlib.sha256(after_quench[0].encode() + after_quench[1].astype('<u4', copy=False).tobytes() + struct.pack('<qqd', *after_quench[2:])).hexdigest()),
     **result_fields(observed_loose))
np.random.set_state(after_quench)

original_random_vector = rotations.vec_random_ndim


def random_vector(shape):
    before = state()
    replay = np.random.RandomState()
    replay.set_state(np.random.get_state())
    raw = replay.normal(size=shape)
    value = original_random_vector(shape)
    emit('global_eigenvector_draw', raw=raw.tolist(), raw_hash=digest(raw),
         normalized=value.tolist(), normalized_hash=digest(value), before_rng=before, after_rng=state())
    return value


rotations.vec_random_ndim = random_vector
original_eig_init = FindLowestEigenVector.__init__
original_eig_run = FindLowestEigenVector.run
original_transverse = FindTransitionState._get_lowest_transverse_eigenvector
phase = 'primary'
eig_count = 0


def eig_init(self, coords, potential, *args, **kwargs):
    global eig_count
    eig_count += 1
    self.audit_id = eig_count
    self.audit_phase = phase
    emit('lowest_eigenvector_input', id=eig_count, phase=phase, coords=coords.tolist(),
         coords_hash=digest(coords), rng=state(),
         dx=float(kwargs.get('dx', 1e-6)), first_order=bool(kwargs.get('first_order', True)))
    original_eig_init(self, coords, potential, *args, **kwargs)
    if self.audit_id == 1 or phase == 'secondary':
        emit('rayleigh_initial', id=self.audit_id, phase=phase,
             coords=self.minimizer.X.tolist(), coords_hash=digest(self.minimizer.X),
             eigenval=float(self.minimizer.energy), rms=float(self.minimizer.rms),
             gradient=self.minimizer.G.tolist(), gradient_hash=digest(self.minimizer.G))

        def accepted_step(**values):
            vector = values['coords']
            emit('rayleigh_step', id=self.audit_id, phase=self.audit_phase,
                 step=int(self.minimizer.iter_number) + 1,
                 eigenval=float(values['energy']), rms=float(values['rms']),
                 coords=vector.tolist(), coords_hash=digest(vector))

        self.minimizer.events.append(accepted_step)


def eig_run(self, niter=None):
    result = original_eig_run(self, niter)
    projection = self.eigpot.orthogZeroEigs
    projected = result.eigenvec.copy()
    if projection is not None:
        projected = projection(projected, self.eigpot.coords)
    leakage = float(np.linalg.norm(result.eigenvec - projected))
    emit('lowest_eigenvector_result', id=self.audit_id, phase=self.audit_phase,
         budget=niter, configured_budget=self.minimizer.nsteps,
         leakage=leakage, **result_fields(result), rng=state())
    return result


def serial(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


helper_source, helper_start = inspect.getsourcelines(original_transverse)
central_line = next(helper_start + index for index, line in enumerate(helper_source)
                    if 'if (not np.all(np.isfinite(restricted))' in line)


def capture_validation(frame, event, arg):
    if frame.f_code != original_transverse.__code__:
        return None
    if event == 'return' or (event == 'line' and frame.f_lineno == central_line):
        keys = ('norm', 'primary_curvature', 'primary_residual', 'curvature',
                'residual', 'coupling', 'restricted', 'direction', 'guess')
        data = {key: serial(frame.f_locals[key]) for key in keys if key in frame.f_locals}
        emit('validation_return' if event == 'return' else 'secondary_central', **data)
    return capture_validation


def transverse(self, coords, gradient):
    global phase
    phase = 'secondary'
    emit('transverse_validation_start', coords=coords.tolist(), coords_hash=digest(coords))
    previous_trace = sys.gettrace()
    sys.settrace(capture_validation)
    try:
        result = original_transverse(self, coords, gradient)
        emit('transverse_validation_result', **result_fields(result))
        return result
    finally:
        sys.settrace(previous_trace)
        phase = 'primary'


FindLowestEigenVector.__init__ = eig_init
FindLowestEigenVector.run = eig_run
FindTransitionState._get_lowest_transverse_eigenvector = transverse


def outer_event(**values):
    coords = values.pop('coords')
    emit('outer_step', **values, coords=coords.tolist(), coords_hash=digest(coords))


wrapped = _PotWrapper(pot)
opt = FindTransitionState(loose.coords, wrapped, event=outer_event)
with open(os.devnull, 'w') as sink, contextlib.redirect_stdout(sink):
    result = opt.run()
emit('transition_state_result', counted_nfev=int(wrapped.nfev), **result_fields(result), rng=state())

# Independent analytic curvature is evaluated only after the algorithm has returned.
energy, gradient, hessian = pot.getEnergyGradientHessian(result.coords)
centered = result.coords.reshape(-1, 3)
centered = centered - centered.mean(axis=0)
rigid = np.column_stack([np.tile(axis, len(centered)) for axis in np.eye(3)] +
                        [np.cross(np.tile(axis, (len(centered), 1)), centered).ravel() for axis in np.eye(3)])
basis = np.linalg.qr(rigid, mode='complete')[0][:, 6:]
values = np.linalg.eigvalsh(basis.T @ hessian @ basis)
emit('final_physics', energy=float(energy), rms=float(np.linalg.norm(gradient) / np.sqrt(gradient.size)),
     rigid_rank=int(np.linalg.matrix_rank(rigid)), physical_lowest=values[:8].tolist(),
     physical_index=int(np.count_nonzero(values < -1e-5)))

# Replay only the returned point, changing one numerical setting at a time.
source = textwrap.dedent(inspect.getsource(original_transverse))
old_project = """    def project(vector, position):
        vector = physical(vector)
        return vector - np.dot(vector, lowest) * lowest
"""
new_project = """    def project(vector, position):
        original = vector
        vector = physical(vector)
        original[:] = vector - np.dot(vector, lowest) * lowest
        return original
"""
assert old_project in source
source = source.replace(old_project, new_project)
source = source.replace('project(product, coords)', 'project(product.copy(), coords)')
namespace = original_transverse.__globals__.copy()
exec(source, namespace)
inplace_helper = namespace['_get_lowest_transverse_eigenvector']
for label, parameters, helper in (
    ('baseline', {}, original_transverse),
    ('budget500', {'nsteps': 500}, original_transverse),
    ('central100', {'first_order': False}, original_transverse),
    ('inplace100', {}, inplace_helper),
    ('central500', {'first_order': False, 'nsteps': 500}, original_transverse),
):
    replay = FindTransitionState(result.coords, pot, eigenvec0=result.eigenvec.copy(),
                                 lowestEigenvectorQuenchParams=parameters, verbosity=0)
    phase = 'secondary:' + label
    fixed = helper(replay, result.coords.copy(), gradient.copy())
    emit('fixed_candidate_replay', label=label, counted_nfev=int(replay.nfev), **result_fields(fixed))
    phase = 'primary'
out.close()
