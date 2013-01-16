import numpy as np
from copy import deepcopy
from pygmin.angleaxis.molecules import create_water
from pygmin.angleaxis import RBSystem
import pygmin.angleaxis.aamindist as am
import time
import gmin_ as GMIN
from pygmin.potentials import GMINPotential

GMIN.initialize()
pot = GMINPotential(GMIN)
nrigid = 10
water = create_water()
system = RBSystem()
system.add_sites([deepcopy(water) for i in xrange(nrigid)])

x = pot.getCoords()
zev = system.zeroEV(x)

eps = 1e-5
for dx in zev:    
    print "ev test", (pot.getEnergy(x) - pot.getEnergy(x + eps*dx))/eps  

dx = np.random.random(x.shape)
dx/=np.linalg.norm(dx)
print "ev test", (pot.getEnergy(x) - pot.getEnergy(x + eps*dx))/eps  


t0 = time.time()
for i in xrange(1000):
    coords1 = np.random.random(6*nrigid)*4
    coords2 = np.random.random(6*nrigid)*4
    
    coords1[:3*nrigid]=0
    coords2[:3*nrigid]=0
    
    measure1 = am.MeasureAngleAxisCluster(system)
    measure2 = am.MeasureRigidBodyCluster(system)
    
    #print 
    measure1.get_dist(coords1, coords2)

t1 = time.time()
for i in xrange(1000):
    coords1 = np.random.random(6*nrigid)*4
    coords2 = np.random.random(6*nrigid)*4
    
    coords1[:3*nrigid]=0
    coords2[:3*nrigid]=0
    
    measure1 = am.MeasureAngleAxisCluster(system)
    measure2 = am.MeasureRigidBodyCluster(system)
    measure2.get_dist(coords1, coords2)
t2 = time.time()
print t1-t0, t2-t1