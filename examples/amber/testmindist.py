from pygmin.storage import Database
from amberSystem import AMBERsystem

natoms = 22
sys = AMBERsystem(natoms)

db = Database("aladipep_database")

m1, m2 = db.minima()[:2]

mindist = sys.get_mindist()
dist, c1, c2 = mindist(m1.coords, m2.coords)
print "distance", dist