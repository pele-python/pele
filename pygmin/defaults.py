from pygmin.optimize import quench
from pygmin.NEB import dimer

quenchRoutine = quench.lbfgs_py
quenchParams = dict()

tsSearchRoutine = dimer
tsSearchParams = dict()

lowestEigenvectorQuenchParams = dict()

NEBquenchRoutine = quench.lbfgs_py
NEBquenchParams = dict()