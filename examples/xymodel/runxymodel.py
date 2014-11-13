import numpy as np
from numpy import cos, sin
from copy import copy
from pele.potentials.xyspin import XYModel


def angle2vec(a):
    return np.array([cos(a), sin(a)])


def printspins(fout, pot, angles):
    for node in pot.G.nodes():
        i = pot.indices[node]
        s = angle2vec(angles[i])
        fout.write("%g %g " % (node[0], node[1]))
        fout.write("%g %g\n" % (s[0], s[1]))


pi = np.pi
L = 24
nspins = L ** 2

pot = XYModel(dim=[L, L], phi=np.pi)

angles = np.random.uniform(-pi, pi, nspins)
print angles

e = pot.getEnergy(angles)
print "energy ", e



# try a quench
if False:
    from pele.optimize import mylbfgs

    ret = mylbfgs(angles, pot)

    print ret


# set up and run basin hopping

from pele.basinhopping import BasinHopping
from pele.takestep.displace import RandomDisplacement
from pele.takestep.adaptive import AdaptiveStepsize
from pele.storage import savenlowest

# should probably use a different take step routine  which takes into account
# the cyclical periodicity of angles
takestep = RandomDisplacement(stepsize=np.pi / 4)
takestepa = AdaptiveStepsize(takestep, frequency=20)
storage = savenlowest.SaveN(500)

bh = BasinHopping(angles, pot, takestepa, temperature=1.01, storage=storage)
bh.run(400)

print "minima found"
with open("out.spin", "w") as fout:
    for min in storage.data:
        print "energy", min.energy
        fout.write("# %g\n" % (min.energy))
        printspins(fout, pot, min.coords)
        fout.write('\n\n')
        """
        view this in gnuplot with the command
        set size ratio -1
        plot 'out.spin' index 0 u 1:2 w p pt 5, '' index 0 u 1:2:($3*0.5):($4*0.5) w vectors
        """

with open("out.energies", "w") as fout:
    for min in storage.data:
        fout.write("%g\n" % (min.energy))

try:
    lowest = storage.data[0].coords
    import matplotlib.pyplot as plt

    x = [node[0] for node in pot.G.nodes()]
    y = [node[1] for node in pot.G.nodes()]
    ilist = [pot.indices[node] for node in pot.G.nodes()]
    v0 = cos(lowest)
    v1 = sin(lowest)
    plt.quiver(x, y, v0, v1)
    a = plt.gca()
    a.set_xlim([-1, max(x) + 1])
    a.set_ylim([-1, max(y) + 1])
    plt.show()
except:
    print "problem ploting with matplotlib"


