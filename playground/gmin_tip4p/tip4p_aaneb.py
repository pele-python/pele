import numpy as np
from pele.potentials import GMINPotential
import gmin_ as GMIN
from pele.utils import rotations
from pele.takestep import buildingblocks
from pele.transition_states import NEB
import pylab as pl
from copy import deepcopy, copy
from pele.optimize import mylbfgs
from pele.angleaxis import RigidFragment, RBSystem
import tip4p
from tip4p import dump_path

# initialize GMIN and potential
GMIN.initialize()
pot = GMINPotential(GMIN)
coords = pot.getCoords()

nrigid = coords.size / 6
print "I have %d water molecules in the system"%nrigid

water = tip4p.water()
#water.S = np.identity(3, dtype="float64")
#water.S*=10.

# define the whole water system
system = RBSystem()
system.add_sites([deepcopy(water) for i in xrange(nrigid)])

# this is an easy access wrapper for coordinates array
ca = system.coords_adapter(coords)

#buildingblocks.rotate(3.0, ca.rotRigid[-1:])
#ret = mylbfgs(coords, pot.getEnergyGradient, iprint=0)
#coords2 = ret[0]
#np.savetxt("coords1_2.txt", coords1)
#np.savetxt("coords2_2.txt", coords2)

coords1 = np.loadtxt("coords1.txt")
coords2 = np.loadtxt("coords2.txt")
print pot.getEnergy(coords1), pot.getEnergy(coords2)

NEBquenchParams = dict()
NEBquenchParams["nsteps"] = 200
NEBquenchParams["iprint"] = 1
NEBquenchParams["maxstep"] = 0.1
NEBquenchParams["maxErise"] = 0.1
NEBquenchParams["tol"] = 1e-6
NEBquenchRoutine = mylbfgs
decp = dict()
decp["local_connect_params"] = dict()
decp["local_connect_params"]["NEBparams"] = dict()
decp["local_connect_params"]["NEBparams"]["NEBquenchParams"] = NEBquenchParams
decp["local_connect_params"]["NEBparams"]["NEBquenchRoutine"] = NEBquenchRoutine

k = 10.
nimages=50
dneb=True

#print coords2[-6:],coords1[-6:]   

path = tip4p.get_path(system, coords1, coords2, nimages)
path_energy = [pot.getEnergy(coords) for coords in path]
# try the old neb
dump_path("interpolate.xyz", system, path)

neb = NEB(path, pot, k=k, dneb=dneb, with_springenergy = False)
neb.optimize()
neb_1 = neb.copy()
dump_path("neb1.xyz", system, neb_1.coords)
neb.optimize()
neb_2 = neb
dump_path("neb2.xyz", system, neb_2.coords)

# try the new neb
aaneb = NEB(path, pot, distance=system.neb_distance, k=k/20., dneb=dneb, with_springenergy = False)
aaneb.optimize()
aaneb_1 = aaneb.copy()
dump_path("aaneb1.xyz", system, aaneb_1.coords)

aaneb.optimize()
aaneb_2 = aaneb
dump_path("aaneb2.xyz", system, aaneb_2.coords)

fig = pl.figure()
nebsteps = NEBquenchParams["nsteps"] 
pl.plot(path_energy, "k", label="interpolation")
pl.plot(neb_1.energies, "b-", label="NEB - %d quenches"%(nebsteps))
pl.plot(neb_2.energies, "b--", label="NEB - %d quenches"%(2*nebsteps))
pl.plot(aaneb_1.energies, "r-", label="AANEB - %d quenches"%(nebsteps))
pl.plot(aaneb_2.energies, "r--", label="AANEB - %d quenches"%(2*nebsteps))
pl.xlabel("neb image")
pl.ylabel("energy")
pl.legend(loc='best')
 
fig2 = pl.figure()

pl.plot(tip4p.path_dist_cart(system, path), "k", label="interpoaltion")
pl.plot(tip4p.path_dist_cart(system, neb_1.coords), "b-", label="NEB 1")
pl.plot(tip4p.path_dist_cart(system, neb_2.coords), "b--", label="NEB 2")
pl.plot(tip4p.path_dist_cart(system, aaneb_1.coords), "r-", label="AANEB 1")
pl.plot(tip4p.path_dist_cart(system, aaneb_2.coords), "r--", label="AANEB 2")

pl.xlabel("neb image")
pl.ylabel("distance")
pl.title("cartesian distances")
pl.legend(loc='best')

fig3 = pl.figure()

pl.plot(tip4p.path_dist_aa(system, path), "k", label="interpoaltion")
pl.plot(tip4p.path_dist_aa(system, neb_1.coords), "b-", label="NEB 1")
pl.plot(tip4p.path_dist_aa(system, neb_2.coords), "b--", label="NEB 2")
pl.plot(tip4p.path_dist_aa(system, aaneb_1.coords), "r-", label="AANEB 1")
pl.plot(tip4p.path_dist_aa(system, aaneb_2.coords), "r--", label="AANEB 2")

pl.xlabel("neb image")
pl.ylabel("distance")
pl.title("aa distances")
pl.legend(loc='best')

pl.show()


