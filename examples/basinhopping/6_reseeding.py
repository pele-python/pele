"""
Example 5: if the energy doesn't improve after 20 basinhopping steps
then do a short monte carlo run at a very high temperature.
"""
from pele.systems import LJCluster
from pele.takestep import RandomDisplacement, Reseeding
from pele.mc import MonteCarlo


class TakeStepMonteCarlo(object):
    def __init__(self, pot, T=10., nsteps=100, stepsize=0.1):
        self.potential = pot
        self.T = T
        self.nsteps = nsteps

        self.mcstep = RandomDisplacement(stepsize=stepsize)

    def takeStep(self, coords, **kwargs):
        # ake a new monte carlo class
        mc = MonteCarlo(coords, self.potential, self.mcstep,
                        temperature=self.T, outstream=None)
        mc.run(self.nsteps)
        coords[:] = mc.coords[:]

    def updateStep(self, acc, **kwargs):
        pass


natoms = 12
niter = 100
system = LJCluster(natoms)

# define a custom takestep routine
potential = system.get_potential()
reseed = TakeStepMonteCarlo(potential, T=100, nsteps=1000)
takestep = system.get_takestep()
stepGroup = Reseeding(takestep, reseed, maxnoimprove=20)

db = system.create_database()
bh = system.get_basinhopping(database=db, takestep=stepGroup)
bh.run(niter)
print "the lowest energy found after", niter, " basinhopping steps is", db.minima()[0].energy
print ""


