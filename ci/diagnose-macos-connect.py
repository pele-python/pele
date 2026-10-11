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
from pele.transition_states import NEBDriver
from pele.landscape.tests.test_double_ended_connect import TestDoubleEndedConnect
from pele.transition_states.tests.test_NEB import _x1, _x2

logging.basicConfig(level=logging.INFO)
trace = open("connect-trace.jsonl", "w", encoding="utf-8")


def record(kind, **data):
    trace.write(json.dumps({"kind": kind, **data}) + "\n")
    trace.flush()
    print(kind, {key: value for key, value in data.items()
                 if key not in ("coords", "path", "energies")}, flush=True)


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


def traced_refine(*args, **kwargs):
    result = refine(*args, **kwargs)
    for label, ret in zip(("ts", "minimum1", "minimum2"), result[1:]):
        if ret is not None:
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
    trace.close()
