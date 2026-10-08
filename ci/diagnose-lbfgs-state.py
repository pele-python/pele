import argparse
import importlib.metadata
import json
import os
import platform
import sys

import numpy as np
import pele
from pele.optimize import LBFGS
from pele.systems import LJCluster

FIELDS = ('s', 'y', 'rho', 'dXold', 'dGold')

# Diagnostic only. Exact continuation equality is checked after all controls.
# Run outside a source checkout with OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=1.

def metric(a, b):
    a, b = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    av, bv = a.view(np.uint64), b.view(np.uint64)
    sign = np.uint64(1 << 63)
    ak = np.where(av & sign, ~av, av | sign)
    bk = np.where(bv & sign, ~bv, bv | sign)
    return {'equal': bool(np.array_equal(a, b)), 'max_abs': float(np.max(np.abs(a - b))),
            'max_ulp': int(np.max(np.maximum(ak, bk) - np.minimum(ak, bk)))}

def align(a):
    return int(a.ctypes.data % 64)

def aligned_copy(a, residue):
    buffer = np.empty(a.size + 8, dtype=np.float64)
    offset = ((residue - buffer.ctypes.data) % 64) // 8
    result = buffer[offset:offset + a.size].reshape(a.shape)
    result[:] = a
    assert align(result) == residue
    return result

def probe(seed, residue=None, cached=False, trace=False):
    np.random.seed(seed)
    system = LJCluster(13)
    potential = system.get_potential()
    initial_coords = system.get_random_configuration()
    original = LBFGS(initial_coords, potential)
    for _ in range(10):
        original.one_iteration()
    state = original.get_state()
    state_values = {field: getattr(state, field).copy() for field in FIELDS}
    checkpoint_coords = original.X.copy()
    checkpoint_gradient = original.G.copy()
    checkpoint_energy = original.energy
    checkpoint_rms = original.rms
    snapshot_alias = {field: bool(np.shares_memory(getattr(original, field), getattr(state, field))) for field in FIELDS}
    checkpoint_alignment = {field: [align(getattr(original, field)), align(getattr(state, field))] for field in FIELDS}
    old_get_step = original._get_LBFGS_step
    original_steps = []
    def original_step(gradient):
        step = old_get_step(gradient)
        original_steps.append((step.copy(), gradient.copy(),
                               {field: align(getattr(original, field)) for field in FIELDS},
                               align(gradient), align(step)))
        return step
    if trace:
        original._get_LBFGS_step = original_step
    for _ in range(10):
        original.one_iteration()
    assert all(np.array_equal(state_values[field], getattr(state, field)) for field in FIELDS)
    kwargs = {'energy': checkpoint_energy, 'gradient': checkpoint_gradient.copy()} if cached else {}
    restarted = LBFGS(checkpoint_coords, potential, **kwargs)
    recomputed = {'energy': metric(checkpoint_energy, restarted.energy),
                  'gradient': metric(checkpoint_gradient, restarted.G),
                  'rms': metric(checkpoint_rms, restarted.rms)}
    if residue is not None:
        state = state._replace(**{field: aligned_copy(getattr(state, field), residue) for field in FIELDS})
        restarted.G = aligned_copy(restarted.G, residue)
    restarted.set_state(state)
    restored = {field: metric(state_values[field], getattr(restarted, field)) for field in FIELDS}
    assert all(value['equal'] for value in restored.values())
    assert restarted.k == state.k and restarted.H0 == state.H0
    assert restarted._have_dXold == state.have_dXold
    restore_alignment = {field: align(getattr(restarted, field)) for field in FIELDS}
    set_alias = {field: bool(np.shares_memory(getattr(restarted, field), getattr(state, field))) for field in FIELDS}
    new_get_step = restarted._get_LBFGS_step
    restart_steps = []
    def restart_step(gradient):
        step = new_get_step(gradient)
        restart_steps.append((step.copy(), gradient.copy(),
                              {field: align(getattr(restarted, field)) for field in FIELDS},
                              align(gradient), align(step)))
        return step
    if trace:
        restarted._get_LBFGS_step = restart_step
    for _ in range(10):
        restarted.one_iteration()
    state1, state2 = original.get_state(), restarted.get_state()
    record = {'seed': seed, 'residue': residue, 'cached': cached,
              'coords': metric(original.X, restarted.X),
              'energy': metric(original.energy, restarted.energy),
              'gradient': metric(original.G, restarted.G), 'rms': metric(original.rms, restarted.rms),
              'history': {field: metric(getattr(state1, field), getattr(state2, field)) for field in FIELDS},
              'H0': metric(state1.H0, state2.H0), 'k': [state1.k, state2.k],
              'have_dXold': [state1.have_dXold, state2.have_dXold],
              'restored_memory': restored,
              'snapshot_alias_to_original': snapshot_alias,
              'set_state_alias_to_snapshot': set_alias,
              'checkpoint_alignment': checkpoint_alignment,
              'restore_alignment': restore_alignment, 'recomputed_checkpoint': recomputed}
    if trace:
        record['steps'] = [{'iteration': i + 1, 'step': metric(a[0], b[0]),
                            'input_gradient': metric(a[1], b[1]),
                            'history_alignment': [a[2], b[2]],
                            'gradient_alignment': [a[3], b[3]],
                            'step_alignment': [a[4], b[4]]}
                           for i, (a, b) in enumerate(zip(original_steps, restart_steps))]
    if not exact(record):
        record['initial_coords'] = initial_coords.tolist()
        record['checkpoint_coords'] = checkpoint_coords.tolist()
        record['final_coords'] = [original.X.tolist(), restarted.X.tolist()]
    return record

def exact(record):
    return (all(record[key]['equal'] for key in ('coords', 'energy', 'gradient', 'rms', 'H0'))
            and all(record['history'][field]['equal'] for field in FIELDS)
            and record['k'][0] == record['k'][1]
            and record['have_dXold'][0] == record['have_dXold'][1])

parser = argparse.ArgumentParser()
parser.add_argument('--seeds', type=int, default=200)
parser.add_argument('--alignments', action='store_true')
parser.add_argument('--output', help='Optional JSON evidence file')
args = parser.parse_args()
environment = {'python': sys.version, 'numpy': np.__version__, 'pele': pele.__file__,
               'platform': platform.platform(),
               'dependencies': {name: importlib.metadata.version(name) for name in ('numpy', 'scipy', 'pele')},
               'threads': {name: os.environ.get(name) for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS')},
               'runtime_paths': {name: os.environ.get(name) for name in ('CONDA_PREFIX', 'LD_LIBRARY_PATH', 'DYLD_LIBRARY_PATH')}}
if sys.platform == "darwin":
    import ctypes
    dyld = ctypes.CDLL(None)
    dyld._dyld_image_count.restype = ctypes.c_uint32
    dyld._dyld_get_image_name.argtypes = [ctypes.c_uint32]
    dyld._dyld_get_image_name.restype = ctypes.c_char_p
    environment["loaded_numeric_libraries"] = [
        dyld._dyld_get_image_name(i).decode()
        for i in range(dyld._dyld_image_count())
        if any(name in dyld._dyld_get_image_name(i).decode().lower()
               for name in ("accelerate", "blas", "lapack", "veclib"))
    ]
print('ENVIRONMENT ' + json.dumps(environment), flush=True)
np.show_config()
records = []
first_mismatch = None
def record_case(record):
    global first_mismatch
    records.append(record)
    if not exact(record) and first_mismatch is None:
        first_mismatch = record
        print('FIRST_MISMATCH ' + json.dumps(record), flush=True)

for seed in range(args.seeds):
    record_case(probe(seed))
if args.alignments:
    for seed in range(min(args.seeds, 10)):
        for residue in range(0, 64, 8):
            record_case(probe(seed, residue=residue, trace=True))
        record_case(probe(seed, cached=True, trace=True))
if first_mismatch is not None:
    seed = first_mismatch['seed']
    for cached in (False, True):
        print('FIRST_SEED_TRACE ' + json.dumps(probe(seed, cached=cached, trace=True)), flush=True)
    for residue in range(0, 64, 8):
        print('FIRST_SEED_ALIGNMENT_TRACE ' + json.dumps(probe(seed, residue=residue, trace=True)), flush=True)
summary = {'cases': len(records), 'exact_mismatches': sum(not exact(r) for r in records),
      'coords_different': sum(not r['coords']['equal'] for r in records),
      'energy_different': sum(not r['energy']['equal'] for r in records),
      'gradient_different': sum(not r['gradient']['equal'] for r in records),
      'checkpoint_gradient_different': sum(not r['recomputed_checkpoint']['gradient']['equal'] for r in records),
      'checkpoint_energy_different': sum(not r['recomputed_checkpoint']['energy']['equal'] for r in records),
      'rms_different': sum(not r['rms']['equal'] for r in records),
      'history_different': sum(any(not r['history'][field]['equal'] for field in FIELDS) for r in records),
      'H0_different': sum(not r['H0']['equal'] for r in records),
      'k_different': sum(r['k'][0] != r['k'][1] for r in records),
      'snapshot_alias_to_original': sum(any(r['snapshot_alias_to_original'].values()) for r in records),
      'set_state_alias_to_snapshot': sum(all(r['set_state_alias_to_snapshot'].values()) for r in records),
      'max_coords_abs': max(r['coords']['max_abs'] for r in records),
      'max_coords_ulp': max(r['coords']['max_ulp'] for r in records)}
print('SUMMARY ' + json.dumps(summary), flush=True)
print('SAMPLE_SEED0 ' + json.dumps(records[0]), flush=True)
if args.output:
    with open(args.output, 'w') as handle:
        json.dump({'environment': environment, 'summary': summary, 'cases': records}, handle, indent=2)
assert first_mismatch is None, 'Exact LBFGS continuation mismatch; see FIRST_MISMATCH diagnostics'
