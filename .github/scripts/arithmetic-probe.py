"""Temporary alignment observer: installed wheel only; no solver changes."""
import contextlib
import hashlib
import io
import json
import math
import os
import platform
import sys
from unittest.mock import patch

import numpy as np
import scipy
from threadpoolctl import threadpool_info

import pele
from pele.potentials import LJ
from pele.transition_states import orthogopt
from pele.transition_states._find_lowest_eig import LowestEigPot

# First start coordinates in the Intel trajectory-1 trace, preserved as binary64.
COORDS = np.array([
    0.4060698174435644, 0.6898899526641312, 0.45686183979755496,
    0.18030109148646645, -0.30599352001143804, 0.6150198737306561,
    -0.21995206491641303, 1.4309439360687817, 1.6964383141250432,
    -0.3699087724010796, 1.3103905196081673, 0.20680306983995178,
    0.45784044301167437, 1.4596644616116117, -1.622801864114092,
    -1.5133405254117345, -1.7594774330948544, 1.2193371355508458,
    1.014335096306328, 1.3562794474775164, 1.7532920395249192,
    1.1218865928542838, -0.14058921192244217, 1.0267541585443247,
    -1.3988685311464435, 0.514255378310923, -1.3081409934979429,
    1.6265929723140875, 0.07943931595738152, -0.3087060913109502,
    -1.1298975531205895, 0.9125380916499378, -0.33969625210588617,
    0.030960445628205567, -1.82929507431884, 0.3475605543088628,
    0.4013239616980934, 0.4347394688059725, 1.621501894198803,
    0.6957665177556281, -0.5330906420321455, -0.27783858193571415,
    0.9423551834691443, -1.5334503730331672, 0.6927290903648889,
    0.6263119739956381, -1.0594941882709423, -1.3546629220206672,
    -0.7508403350744097, -0.5190047589775054, 0.22830553758379873,
    -0.43551805621540174, 1.8902502289161742, -1.407366978117771,
])
# Seed 0, consume LJ18's 54 uniform draws, then the first 54 Gaussian draws.
GAUSSIAN = np.array([
    -0.5096521817516535, -0.4380743016111864, -1.2527953600499262,
    0.7774903558319101, -1.6138978475579515, -0.2127402802139687,
    -0.8954665611936756, 0.386902497859262, -0.510805137568873,
    -1.180632184122412, -0.028182228338654868, 0.42833187053041766,
    0.06651722238316789, 0.3024718977397814, -0.6343220936809636,
    -0.3627411659871381, -0.672460447775951, -0.3595531615405413,
    -0.813146282044454, -1.7262826023316769, 0.17742614225375283,
    -0.4017809362082619, -1.6301983469660446, 0.4627822555257742,
    -0.9072983643832422, 0.05194539579613895, 0.7290905621775369,
    0.12898291075741067, 1.1394006845433007, -1.2348258203536526,
    0.402341641177549, -0.6848100909403132, -0.8707971491818818,
    -0.5788496647644155, -0.31155253212737266, 0.05616534222974544,
    -1.1651498407833565, 0.9008264869541871, 0.46566243973045984,
    -1.5362436862772237, 1.4882521937955997, 1.8958891760305832,
    1.1787795711596507, -0.17992483581235091, -1.0707526215105425,
    1.0544517269311366, -0.40317694697317963, 1.2224450703824274,
    0.2082749780768603, 0.9766390364837128, 0.3563663971744019,
    0.7065731681919482, 0.010500020720820478, 1.7858704939058352,
])
assert COORDS.shape == GAUSSIAN.shape == (54,)
output = open(sys.argv[1], "w")


def record(kind, **values):
    output.write(json.dumps(dict(kind=kind, **values)) + "\n")
    output.flush()


def aligned(values, offset):
    storage = np.empty(values.size + 15)
    start = (-storage.ctypes.data // 8) % 8 + offset
    view = storage[start:start + values.size]
    view[:] = values
    assert view.ctypes.data % 64 == offset * 8
    assert view.tobytes() == values.tobytes()
    return view


baseline = {}
first = {}
case = None


def stage(path, name, value):
    array = np.asarray(value)
    assert np.all(np.isfinite(array))
    result = [float(v).hex() for v in array.ravel()]
    key = (path, name)
    reference = baseline.setdefault(key, result)
    changed = result != reference
    witness = None
    if changed:
        index = next(i for i, (a, b) in enumerate(zip(reference, result)) if a != b)
        witness = dict(index=index, baseline=reference[index], observed=result[index])
        first.setdefault(path, dict(stage=name, **witness))
    record("stage", case=case, path=path, stage=name, hex=result,
           address=array.ctypes.data, mod64=array.ctypes.data % 64,
           differs=changed, witness=witness)


runtime = io.StringIO()
with contextlib.redirect_stdout(runtime):
    np.show_config()
    np.show_runtime()
record("runtime", python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
       platform=platform.platform(), pele=pele.__file__, numpy_file=np.__file__,
       numpy_config_runtime=runtime.getvalue(), threadpools=threadpool_info(),
       threads={k: os.environ.get(k) for k in (
           "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
           "MKL_NUM_THREADS", "BLIS_NUM_THREADS")})
record("inputs", coords_hex=[float(v).hex() for v in COORDS],
       gaussian_hex=[float(v).hex() for v in GAUSSIAN],
       coords_sha256=hashlib.sha256(COORDS.tobytes()).hexdigest(),
       gaussian_sha256=hashlib.sha256(GAUSSIAN.tobytes()).hexdigest(),
       fsum_square=math.fsum(float(v) * float(v) for v in GAUSSIAN).hex(),
       fsum_norm=math.sqrt(math.fsum(float(v) * float(v) for v in GAUSSIAN)).hex(),
       fsum_dot_coords=math.fsum(float(v) * float(x) for v, x in zip(GAUSSIAN, COORDS)).hex())
lj = LJ()
_, true_gradient = lj.getEnergyGradient(COORDS)
original_norm = np.linalg.norm

for coord_offset in range(8):
    for vector_offset in range(8):
        case = [coord_offset, vector_offset]
        first = {}
        x = aligned(COORDS, coord_offset)
        v = aligned(GAUSSIAN, vector_offset)
        record("case", case=case, coords_address=x.ctypes.data, vector_address=v.ctypes.data,
               coords_mod64=x.ctypes.data % 64, vector_mod64=v.ctypes.data % 64,
               coords_byte_equal=x.tobytes() == COORDS.tobytes(),
               vector_byte_equal=v.tobytes() == GAUSSIAN.tobytes())
        stage("numpy", "norm", original_norm(v))
        stage("numpy", "dot_self", np.dot(v, v))
        stage("numpy", "dot_coords", np.dot(v, x))
        v /= original_norm(v)
        stage("numpy", "normalized", v)
        for norm in (False, True):
            v = aligned(GAUSSIAN, vector_offset)
            stage("orthogopt", "norm=" + str(norm), orthogopt(v, x, norm=norm))
        energy, gradient = lj.getEnergyGradient(x)
        stage("LJ", "energy", energy)
        stage("LJ", "gradient", gradient)

        for method in ("getEnergy", "getEnergyGradient"):
            path = "LowestEigPot." + method
            norm_count = orthog_count = 0

            def traced_norm(value, *args, **kwargs):
                global norm_count
                norm_count += 1
                stage(path, f"norm{norm_count}.input", value)
                result = original_norm(value, *args, **kwargs)
                stage(path, f"norm{norm_count}.result", result)
                return result

            def traced_orthog(value, coords):
                global orthog_count
                orthog_count += 1
                stage(path, f"orthog{orthog_count}.input", value)
                result = orthogopt(value, coords)
                stage(path, f"orthog{orthog_count}.result", result)
                return result

            class TracedLJ:
                def getEnergyGradient(self, coords):
                    stage(path, "displaced_coords", coords)
                    result = lj.getEnergyGradient(coords)
                    stage(path, "displaced_energy", result[0])
                    stage(path, "displaced_gradient", result[1])
                    return result

            eigpot = LowestEigPot(x, TracedLJ(), orthogZeroEigs=traced_orthog,
                                  gradient=true_gradient, first_order=True)
            # Retain the requested address alignment after the constructor's copy.
            eigpot.coords = x
            v = aligned(GAUSSIAN, vector_offset)
            with patch.object(np.linalg, "norm", traced_norm):
                result = getattr(eigpot, method)(v)
            if method == "getEnergy":
                stage(path, "curvature", result)
            else:
                stage(path, "curvature", result[0])
                stage(path, "gradient", result[1])
        assert x.tobytes() == COORDS.tobytes()
        record("first_difference", case=case, paths=first)
record("complete", cases=64)
output.close()
