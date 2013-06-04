import random
import pickle
from pele.storage.database import Database
#from math import pi
from pele.systems.oxdna import OXDNATakestep, export_xyz, OXDNAScrewStep

# number of trial configurations to try
nconf = 1000
# generate these from the n lowest minima
from_nlowest = 100
# open the database with minima
db = Database(db = "storage.sqlite")

minima = db.minima()

# make sure from_nlowest is not larger than # of minima
from_nlowest = max(from_nlowest, len(minima))

# you can try a different step routine
step = OXDNATakestep(displace=0.0, rotate=0.8, rotate_around_backbone=False)

trial_configurations = []
for i in xrange(nconf):
    x = random.choice(minima[0:nconf])
    coords = x.coords.copy()
    step.takeStep(coords)
    trial_configurations.append(coords)
    
pickle.dump(trial_configurations, open("quench_benchmark.dat", "w"))

fl = open("quench_benchmark.xyz", "w")
for conf in trial_configurations:
    export_xyz(fl, conf)

fl.close()
