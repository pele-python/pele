import numpy as np
from copy import copy
from pygmin.potentials.heisenberg_spin_RA import HeisenbergModelRA
import pygmin.utils.rotations as rotations
from numpy import sin, cos

import networkx as nx

from pygmin.potentials.heisenberg_spin import make3dVector,  make2dVector, coords2ToCoords3, coords3ToCoords2, grad3ToGrad2

def getm(coords2):
    coords3 = coords2ToCoords3( coords2 )
    m = np.linalg.norm( coords3.sum(0) ) / nspins
    return m


def printspins(fout, pot, coords2):
    s = coords2ToCoords3( coords2 )
    h = pot.fields
    c = coords2ToCoords3( coordsinit )
    for node in pot.G.nodes():
        i = pot.indices[node]
        fout.write( "%g %g %g %g %g %g %g %g\n" % (node[0], node[1], \
            s[i,0], s[i,1], s[i,2], h[i,0], h[i,1], h[i,2] ) )




pi = np.pi
L = 4
nspins = L**2

#phases = np.zeros(nspins)
pot = HeisenbergModelRA( dim = [L,L], field_disorder = 2. ) #, phases=phases)

coords = np.zeros([nspins, 2])
for i in range(nspins): 
    vec = rotations.vec_random()
    coords[i,:] = make2dVector(vec)
coords = np.reshape(coords, [nspins*2])
#print np.shape(coords)
coordsinit = np.copy(coords)

#print "fields", pot.fields
print coords


e = pot.getEnergy(coords)
print "energy ", e

print "try a quench"
from pygmin.optimize._quench_new import mylbfgs
ret = mylbfgs(coords, pot)

print "quenched e = ", ret.energy, "funcalls", ret.nfev
print ret.coords

m = getm( ret[0] )
print "magnetization after quench", m


#do basin hopping

from pygmin.basinhopping import BasinHopping
from pygmin.takestep.displace import RandomDisplacement
from pygmin.takestep.adaptive import AdaptiveStepsize
from pygmin.storage import savenlowest

takestep = RandomDisplacement(stepsize = np.pi/4)
takestepa = AdaptiveStepsize(takestep, frequency = 10)
storage = savenlowest.SaveN(20)

bh = BasinHopping( coords, pot, takestepa, temperature = 1.01, storage = storage)
bh.run(200)

print "lowest structures fount:"
with open("out.spins", "w") as fout:
    for min in storage.data:
        m = getm( min.coords )
        print "energy", min.energy, "magnetization", m
        fout.write( "energy %g magnetization %g\n" % (min.energy, m) )
        printspins(fout, pot, min.coords)
        fout.write("\n\n")

"""
view the spins with gnuplot using the command
h = 2.
s = 0.7
splot 'out.spins' u 1:2:(0) w p pt 5, '' index 1 u 1:2:(0):($6/h):($7/h):($8/h) w vectors t "fields", '' index 1 u 1:2:(0):($3*s):($4*s):($5*s) w vectors t "spins"
"""


