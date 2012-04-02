import numpy as np
import potentials.lj as lj
import mingmin
import saveit
import mingmin as gmin

potential = lj.LJ(1.0, 1.0, None)

coords=np.loadtxt('coords')
# avoid this is this really needed?
natoms=coords.shape[0]
coords = coords.reshape(coords.size)

# run monte carlo
savelowest = saveit.saveit()
# potential, coords, natoms, mcsteps, temperature, stepsize, equi steps, savelowest
gmin.monteCarlo(potential, coords, natoms, 1000, 1.0, 0.3, 100, savelowest)

# print results
with open("lowest", "w") as fout:
  gmin.printxyz(fout,  savelowest.lowestcoords, natoms, savelowest.lowestE)

