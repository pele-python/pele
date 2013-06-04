from pele.systems import LJCluster
import numpy as np
import scipy.sparse
import scikits.sparse.cholmod as cholmod
import time
import pele.transition_states as ts
from pele.utils.hessian import get_sorted_eig

natoms = 1000
system = LJCluster(natoms)
#system.params.structural_quench_params.debug = True
#system.params.structural_quench_params.iprint = 100
db = system.create_database()
bh = system.get_basinhopping(db)
bh.run(1)

m = db.minima()[0]
coords = m.coords

potential = system.get_potential()
energy, gradient, hessian = potential.getEnergyGradientHessian(coords)

dummy_vec = ts.gramm_schmidt(ts.zeroEV_cluster(coords))

shifted_hess = hessian.copy()

for i in range(6):
    shifted_hess += np.outer(dummy_vec[i], dummy_vec[i])

shifted_eval, shifted_evec = get_sorted_eig(shifted_hess)

print "First log sum: ", np.sum(np.log(shifted_eval[6:]))

sparse_hess = scipy.sparse.csc_matrix(shifted_hess)
factor = cholmod.cholesky(sparse_hess)

diagonal = np.diagonal(factor.L().todense())

logar  = 2 * np.log(diagonal)
log_sum = np.sum(logar)
print "Second log sum: ", log_sum