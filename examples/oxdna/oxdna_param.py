
from pygmin.storage.database import Database
import numpy as np
from pygmin.potentials import GMINPotential
import oxdnagmin_ as GMIN
from pygmin import takestep
from math import pi
from pygmin.utils.rbtools import CoordsAdapter
from pygmin.utils import rotations
from pygmin.basinhopping import BasinHopping
import parameters
import time

t0 = time.clock()

# This is the takestep routine for OXDNA. It is a standard rigid body takestep
# routine, but I put it here to be able to start modifying it
class OXDNATakestep(takestep.TakestepInterface):
    def __init__(self, displace=1.0, rotate=0.5*pi):
        self.displace = displace
        self.rotate = rotate
        
    def takeStep(self, coords, **kwargs):
        # easy access to coordinates
        ca = CoordsAdapter(nrigid=coords.size/6, coords = coords)
        
        # random displacement for positions
        ca.posRigid[:] += 2.*self.displace*(np.random.random(ca.posRigid.shape)-0.5)
        
        # random rotation for angle-axis vectors
        takestep.rotate(self.rotate, ca.rotRigid)
        
    # this is necessary for adaptive step taking
    def scale(self, factor):
        self.rotate *= factor
        self.displace *= factor

# this class should generate a fully random configuration
class OXDNAReseed(takestep.TakestepInterface):
    def __init__(self, radius=3.0):
        self.radius = radius
    
    def takeStep(self, coords, **kwargs):
        # easy access to coordinates
        ca = CoordsAdapter(nrigid=coords.size/6, coords = coords)
        
        # random displacement for positions
        #ca.posRigid[:] = 2.*self.radius*(np.random.random(ca.posRigid.shape)-0.5)
        
        # random rotation for angle-axis vectors
        for rot in ca.rotRigid:
            rot[:] = rotations.random_aa()

    def check_converged(E, coords):
        if(E<(parameters.TARGET+parameters.EDIFF)):
                  fl = open("stat.dat", "a")
                  print "#found minimum"
                  t1= time.clock()
                  timespent= t1 - t0
                  fl.write("#quenches, functioncalls, time\n")
                  fl.write("%d %d %f\n"%(opt.stepnum   , potential.ncalls   , timespent))
                  fl.close()
                  exit()

           
# initialize GMIN
GMIN.initialize()
# create a potential which calls GMIN
potential = GMINPotential(GMIN)
# get the initial coorinates
coords=potential.getCoords()
coords=np.random.random(coords.shape)
# create takestep routine

# we combine a normal step taking
step = OXDNATakestep(displace=parameters.displace, rotate=parameters.rotate)
# with a generate random configuration
genrandom = OXDNAReseed()
# in a reseeding takestep procedure
reseed = takestep.Reseeding(step, genrandom, maxnoimprove=parameters.reseed)
    
# store all minima in a database
db = Database(db="storage.sqlite", accuracy=1e-2)

# create Basinhopping object
opt = BasinHopping(coords, potential, reseed, db.minimum_adder())

# run for 100 steps
opt.run(parameters.nsteps)

# now dump all the minima
i=0
for m in db.minima():
    i+=1
    GMIN.userpot_dump("lowest_%03d.dat"%(i), m.coords)
