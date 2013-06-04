import numpy as np
import pickle
from pele.potentials.lj import LJ
from pele.NEB.NEB import NEB
import pylab as pl

dataset = pickle.load(open("coords.3.dat", "r"))

pot = LJ()

for coords1,coords2 in dataset:   
    neb = NEB(coords1,coords2,pot)
    neb.optimize()
    pl.plot(neb.energies)
    pl.show()
