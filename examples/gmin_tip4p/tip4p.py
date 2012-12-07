from math import sin, cos, pi
import numpy as np
from pygmin.potentials import GMINPotential
from pygmin.angleaxis import RigidFragment, RBSystem
import gmin_ as GMIN
from pygmin.utils import rotations
from pygmin.optimize import quench
from pygmin.takestep import buildingblocks
from pygmin.transition_states import NEB, InterpolatedPath
from pygmin import defaults
import pylab as pl
from copy import deepcopy
from pygmin.optimize import quench
from pygmin.mindist.rmsfit import findrotation_kabsch

def dump_path(filename, system, path):
    fl = open(filename, "w")
    lbls = system.get_atom_labels()
    for c in path:
        from pygmin.printing import printAtomsXYZ
        atomistic = system.to_atomistic(c)
        fl.write("%d\n\n"%len(lbls))
        for lbl, x in zip(lbls, atomistic):
            fl.write("%s %f %f %f\n"%(lbl, x[0], x[1], x[2]))
    fl.close()
    
# initialize GMIN and potential
GMIN.initialize()
pot = GMINPotential(GMIN)
coords = pot.getCoords()

nrigid = coords.size / 6
print "I have %d water molecules in the system"%nrigid


# define a water molecule
water = RigidFragment()
rho   = 0.9572
theta = 104.52/180.0*pi      
water.add_atom("O", np.array([0., 0., 0.]), 16.)
water.add_atom("H", rho*np.array([0.0, sin(0.5*theta), cos(0.5*theta)]), 1.)
water.add_atom("H", rho*np.array([0.0, -sin(0.5*theta), cos(0.5*theta)]), 1.)
water.finalize_setup()

#water.S = np.identity(3, dtype="float64")
#water.S*=10.
# define the whole water system
system = RBSystem()
system.add_sites([deepcopy(water) for i in xrange(nrigid)])

# this is an easy access wrapper for coordinates array
ca = system.coords_adapter(coords)

d = []
d2 = []
s = []
p=ca.rotRigid[3].copy()
p/=np.linalg.norm(p)
p*=pi
coords1 = coords.copy()
ca.rotRigid[3]=ca.rotRigid[3] + p
coords2 = coords.copy()

# try rotation symmetry
#E = pot.getEnergy(coords)
#ca.rotRigid[1] = rotations.rotate_aa(np.array([0., 0., pi]), ca.rotRigid[1])
#print "Energy before and after rotation molecule 2 by Pi around its symmetry axis", E, pot.getEnergy(coords)
#
buildingblocks.rotate(3.0, ca.rotRigid)
ret = quench.mylbfgs(coords, pot.getEnergyGradient, iprint=0)
#print ret[1]
coords2 = ret[0]

#np.savetxt("coords1.txt", coords1)
#np.savetxt("coords2.txt", coords2)
coords1 = np.loadtxt("coords1.txt")
coords2 = np.loadtxt("coords2.txt")
defaults.NEBquenchParams["nsteps"] = 500
defaults.NEBquenchParams["iprint"] = -1
defaults.NEBquenchParams["maxstep"] = 0.1
#defaults.NEBquenchParams["maxErise"] = 0.1
defaults.NEBquenchRoutine = quench.fire
k = 100.
dneb=True


c1 = system.coords_adapter(coords1)
c2 = system.coords_adapter(coords2)
R = findrotation_kabsch(c2.posRigid, c1.posRigid)
#R = rotations.aa2mx(p)
for x, p in zip(c2.posRigid, c2.rotRigid):
    x[:] = np.dot(R, x)
    p[:] = rotations.rotate_aa(p, rotations.mx2aa(R))

print coords2[-6:] ,coords1[-6:]   
for p1, p2 in zip(c1.rotRigid,c2.rotRigid):
    theta1 = np.linalg.norm(rotations.rotate_aa(p2,-p1))
    p2n = rotations.rotate_aa(np.array([0., 0., pi]), p2)
    theta2 = np.linalg.norm(rotations.rotate_aa(p2n,-p1))
    theta1 -= int(theta1/2./pi)*2.*pi
    theta2 -= int(theta2/2./pi)*2.*pi
    if(theta2 < theta1): 
        p2[:]=p2n
    print theta2, theta1
#print coords2[-6:],coords1[-6:]   

path=[x for x in InterpolatedPath(coords2, coords1, 30, interpolator=system.interpolate)]
for i in path:
    print pot.getEnergy(i)
print "foo"
system.align_path(path)
print "foo"
dump_path("interpolate.xyz", system, path)

neb = NEB(path, pot, k=k, dneb=dneb, with_springenergy = False)

neb.optimize()
neb.optimize()
dump_path("neb.xyz", system, neb.coords)

neb2 = NEB(path, pot, distance=system.neb_distance, k=k/20., dneb=dneb, with_springenergy = False)
neb2.optimize()
dump_path("neb2_1.xyz", system, neb2.coords)
e21 = neb2.energies.copy()
dist21_cart = []
dist21_aa = []
for i in xrange(0,neb2.nimages-1):
    dist21_cart.append(np.sqrt(np.linalg.norm(neb2.coords[i] - neb2.coords[i+1])**2))
    dist21_aa.append(np.sqrt(system.distance_squared(neb2.coords[i], neb2.coords[i+1])))


neb2.optimize()
dump_path("neb2.xyz", system, neb2.coords)
#print path

import pickle
pickle.dump(neb.coords, open("path_cart.dat", "w"))
pickle.dump(neb2.coords, open("path_aa.dat", "w"))

import pylab as pl
fig = pl.figure()
nebsteps = defaults.NEBquenchParams["nsteps"] 
pl.plot(neb.energies, label="NEB - %d quenches"%(nebsteps))
pl.plot(neb2.energies, label="AANEB - %d quenches"%(nebsteps))
pl.plot(e21, label="AANEB - %d quenches"%(2*nebsteps))
pl.xlabel("neb image")
pl.ylabel("energy")
pl.legend(loc='best')

dist1_cart = []
dist1_aa = []
for i in xrange(0,neb.nimages-1):
    dist1_cart.append(np.sqrt(neb.distance(neb.coords[i], neb.coords[i+1])[0]))
    dist1_aa.append(np.sqrt(system.distance_squared(neb.coords[i], neb.coords[i+1])))

dist2_cart = []
dist2_aa = []
for i in xrange(0,neb2.nimages-1):
    dist2_cart.append(np.sqrt(np.linalg.norm(neb2.coords[i] - neb2.coords[i+1])**2))
    dist2_aa.append(np.sqrt(system.distance_squared(neb2.coords[i], neb2.coords[i+1])))

dist0_cart = []
dist0_aa = []
for i in xrange(0,neb2.nimages-1):
    dist0_cart.append(np.linalg.norm(path[i] - path[i+1])**2)
    dist0_aa.append(system.distance_squared(path[i], path[i+1]))
    
fig2 = pl.figure()

pl.plot(dist0_cart, label="interpoaltion - cartesian distances")
pl.plot(dist1_cart, label="Standard NEB - cartesian distances")
pl.plot(dist2_cart, "--", label="AngleAxis NEB - cartesian distances")
pl.xlabel("neb image")
pl.ylabel("energy")
pl.legend(loc='best')
fig3 = pl.figure()

pl.plot(dist0_aa, label="interpoaltion - aa distances")
pl.plot(dist1_aa, label="Standard NEB - aa distances")
pl.plot(dist21_aa, "--", label="AngleAxis NEB - aa distances 1")
pl.plot(dist2_aa, "--", label="AngleAxis NEB - aa distances")
pl.xlabel("neb image")
pl.ylabel("energy")
pl.legend(loc='best')

pl.show()


