import numpy as np
from pap_example import PAPSystem
import numpy as np


def add(db, system, x):
    pot = system.get_potential()
    x = x.flatten()
    E=pot.getEnergy(x)
    print "minimum before quench %e"%E
    opt = system.get_minimizer()
    ret = opt(x)
    E = ret[1]
    x = ret[0]
    m=db.addMinimum(E, x)
    print "added minimum %d with energy %e"%(m._id, m.energy)
    print

system = PAPSystem()
db = system.create_database("pap.sqlite")

add(db, system, np.loadtxt("coords1"))
add(db, system, np.loadtxt("coords2"))

