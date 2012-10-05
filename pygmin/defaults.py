from optimize import quench
from NEB import dimer

quenchRoutine = quench.lbfgs_py
quenchParams = dict()

tsSearchRoutine = dimer
tsSearchParams = dict()

NEBquenchRoutine = quench.lbfgs_py
NEBquenchParams = dict()