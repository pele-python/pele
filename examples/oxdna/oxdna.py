import numpy as np
import oxdnagmin_ as GMIN
from pygmin.potentials.gminpotential import GMINPotential
import pygmin.basinhopping as bh
from pygmin.takestep import displace
from pygmin.storage.database import Database


# initialize GMIN
GMIN.initialize()
# create a potential which calls GMIN
potential = GMINPotential(GMIN)
# get the initial coorinates
coords=potential.getCoords()

# create takestep routine
step = displace.RandomDisplacement(stepsize=0.5)

# store all minima in a database
db = Database(db="storage.sqlite")

# create Basinhopping object
opt = bh.BasinHopping(coords, potential, step, storage=db.minimum_adder())

# run for 100 steps
opt.run(100)

# now dump all the minima
i=0
for m in db.minima():
    i+=1
    print m.coords
    GMIN.userpot_dump("lowest_%03d.dat"%(i), m.coords)