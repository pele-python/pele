from pele.thermodynamics import normalmode_frequencies, logproduct_freq2
from tip4p_system import TIP4PSystem
from pele.angleaxis.aamindist import TransformAngleAxisCluster
from pele.utils import rotations
import numpy as np
from pele.mindist import PointGroupOrderCluster
from pele.angleaxis import ExactMatchAACluster

system = TIP4PSystem()
db = system.create_database(db="tip4p_8.sqlite")
pot = system.get_potential()
transform = TransformAngleAxisCluster(system.aasystem)
coords = db.minima()[0].coords
#coords = db.transition_states()[0].coords
energy = pot.getEnergy(coords)
#coords = np.loadtxt("2.txt").flatten()
print energy

match = ExactMatchAACluster(system.aasystem)
get_pgorder = PointGroupOrderCluster(match)

coords = db.minima()[0].coords
hess = pot.NumericalHessian(coords, eps=1e-4)
metric = system.aasystem.metric_tensor(coords)

print get_pgorder(coords)
frq = normalmode_frequencies(hess, metric)
print frq
fvib = logproduct_freq2(frq, 6)[1]
print fvib

beta = 1./2.479
for m in db.minima():
    hess = pot.NumericalHessian(m.coords, eps=1e-5)
    metric = system.aasystem.metric_tensor(m.coords)

    pgorder =  get_pgorder(m.coords)
    frq = normalmode_frequencies(hess, metric, eps=1e-3)
    fvib = logproduct_freq2(frq, 6, eps=1e-3)[1]
    print fvib, 0.5*fvib /beta 
    m.pgorder = pgorder
    m.fvib = fvib 

for ts in db.transition_states():
    hess = pot.NumericalHessian(ts.coords, eps=1e-5)
    metric = system.aasystem.metric_tensor(ts.coords)

    pgorder = get_pgorder(ts.coords)
    frq = normalmode_frequencies(hess, metric, eps=1e-3)
    fvib = logproduct_freq2(frq, 6, nnegative=1, eps=1e-3)[1]
    print fvib, 0.5*fvib /beta
    ts.pgorder = pgorder
    ts.fvib = fvib 

db.session.commit()