import numpy as np
from pap_example import PAPSystem

system = PAPSystem()

x1 = np.loadtxt("coords1").flatten()
db = system.create_database("pap.sqlite")
pot = system.get_potential()

E=pot.getEnergy(x1)
db.addMinimum(E, x1)
x1 = np.loadtxt("coords2").flatten()
E=pot.getEnergy(x1)
db.addMinimum(E, x1)
