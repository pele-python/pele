import scipy
import numpy as np
from pele.systems import LJCluster
from pele.mindist import ExactMatchAtomicCluster
from pele.mindist import PointGroupOrderCluster
from pele.thermodynamics import normalmode_frequencies, logproduct_freq2
    
beta = 1.
system = LJCluster(13)

db = system.create_database()
pot = system.get_potential()

bh = system.get_basinhopping(database=db)
bh.run(50)


min1 = db.minima()[0]
coords = min1.coords
print
print "Done with basinghopping, performing frequency analysis"
print
#min1 = db.transition_states()[0]


# determine point group order of system
determine_pgorder = PointGroupOrderCluster(system.get_compare_exact())
pgorder = determine_pgorder(min1.coords)
# free energy from symmetry
Fpg = np.log(pgorder)/beta

# get the hession
e, g, hess = pot.getEnergyGradientHessian(min1.coords)
# TODO: go to reduced coordinates here 

# get the eigenvalues
freqs2 = normalmode_frequencies(hess)

# analyze eigenvalues
n, lnf = logproduct_freq2(freqs2, 6)

Ffrq = (n*np.log(beta) + 0.5*lnf /beta) 

# print a short summary
E = pot.getEnergy(coords)
print "point group order:", pgorder 
print "number of zero eigenvalues:", (np.abs(freqs2) < 1e-4).sum()
print "number of negative eigenvalues:", (freqs2 < -1e-4).sum()
print "free energy at beta=%f:"%beta, E + Ffrq + Fpg
print "contributions from"
print "potential energy", E
print "frequencies", Ffrq
print "point group order", Fpg

