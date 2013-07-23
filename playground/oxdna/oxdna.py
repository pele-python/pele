import numpy as np
import oxdnagmin_ as GMIN
from pele.potentials.gminpotential import GMINPotential
import pele.basinhopping as bh
from pele.takestep import displace
from pele.storage.database import Database


# initialize GMIN
GMIN.initialize()
# create a potential which calls GMIN
potential = GMINPotential(GMIN)
# get the initial coorinates
coords=potential.getCoords()
coords=np.random.random(coords.shape)
# create takestep routine
step = displace.RandomDisplacement(stepsize=1.)

# store all minima in a database
db = Database(db="storage.sqlite", accuracy=1e-2)

# create Basinhopping object
opt = bh.BasinHopping(coords, potential, step, db.minimum_adder())

# run for 100 steps
opt.run(1000)

# now dump all the minima
i=0
for m in db.minima():
    i+=1
    GMIN.userpot_dump("lowest_%03d.dat"%(i), m.coords)