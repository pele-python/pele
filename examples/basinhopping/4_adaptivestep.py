"""
Example 4: modify the parameters of the adaptive stepsize
note that the system class uses adaptive stepsize by default, 
so the previous examples also use adaptive stepsize
"""
from pele.systems import LJCluster
from pele.takestep import RandomDisplacement, AdaptiveStepsizeTemperature

natoms = 12
niter = 100
system = LJCluster(natoms)

db = system.create_database()

# create takestep routine manually
step = RandomDisplacement(stepsize=0.3)
# wrap it in the class which will adaptively improve the stepsize
wrapped_step = AdaptiveStepsizeTemperature(step, interval=30)

bh = system.get_basinhopping(database=db, takestep=wrapped_step)
bh.run(niter)
print "the lowest energy found after", niter, " basinhopping steps is", db.minima()[0].energy
print ""
