"""
Example 5: Adding a custom takestep routine.  This example
takes 100 monte carlo steps as one basin hopping step
"""
from pele.systems import LJCluster
from pele.takestep import RandomDisplacement
from pele.mc import MonteCarlo


class TakeStepMonteCarlo(object):
    def __init__(self, pot, T=10., nsteps=100, stepsize=0.1):
        self.potential = pot
        self.T = T
        self.nsteps = nsteps

        self.mcstep = RandomDisplacement(stepsize=stepsize)

    def takeStep(self, coords, **kwargs):
        # make a new monte carlo class
        mc = MonteCarlo(coords, self.potential, self.mcstep,
                        temperature=self.T, outstream=None)
        mc.run(self.nsteps)
        coords[:] = mc.coords[:]

    def updateStep(self, acc, **kwargs):
        pass


natoms = 12
niter = 100
system = LJCluster(natoms)

db = system.create_database()

# create takestep routine manually
potential = system.get_potential()
step = TakeStepMonteCarlo(potential)

bh = system.get_basinhopping(database=db, takestep=step)
bh.run(niter)
print "the lowest energy found after", niter, " basinhopping steps is", db.minima()[0].energy
print ""
