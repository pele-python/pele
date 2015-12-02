import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.misc import comb, factorial
from pele.optimize._quench import lbfgs_cpp, modifiedfire_cpp
from pele.potentials import MeanFieldPSpinSpherical
from joblib import Parallel, delayed


p=3
nspins=20
interactions = np.ones(np.power(nspins,p))
coords = np.ones(nspins)
pot = MeanFieldPSpinSpherical(interactions, nspins, p, tol=1e-6)
e = pot.getEnergy(coords)
assert e == -comb(nspins,p)/nspins

interactions = np.random.normal(0, np.sqrt(factorial(p)), [nspins for i in xrange(p)])
#interactions =  (interactions + interactions.transpose())/2
for i in xrange(nspins):
    for j in xrange(i, nspins):
        for k in xrange(j, nspins):
            interactions[k][i][j] = interactions[i][j][k]
            interactions[k][j][i] = interactions[i][j][k]
            interactions[j][k][i] = interactions[i][j][k]
            interactions[i][k][j] = interactions[i][j][k]
            interactions[j][i][k] = interactions[i][j][k]
interactions = interactions.flatten()
pot = MeanFieldPSpinSpherical(interactions, nspins, p, tol=1e-6)
energies = []
nfevs = []
coords_list = []

for _ in xrange(1000):
    coords = np.random.normal(0, 1, nspins)
    coords /= (np.linalg.norm(coords)/np.sqrt(nspins))
    coords_list.append(coords)

def func(coords):
    #print coords
    print "start energy", pot.getEnergy(coords)
    results = lbfgs_cpp(coords, pot, nsteps=1e5, tol=1e-5, iprint=-1, maxstep=10)                                                                                                                                 
    #results = modifiedfire_cpp(coords, pot, nsteps=1e6, tol=1e-5, iprint=-1)
    print "quenched energy", results.energy
    if results.success:
        return [results.energy, results.nfev]

out = Parallel(n_jobs=max(1,7))(delayed(func)(x) for x in coords_list)
out = np.array(out)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(out[:,0]/nspins)
ax.set_xlabel('E/N')
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.hist(out[:,1])
ax1.set_xlabel('nfev')
fig.savefig('energy_histogram_n{}.pdf'.format(nspins))
fig1.savefig('nfev_histogram_n{}.pdf'.format(nspins))
plt.show()
