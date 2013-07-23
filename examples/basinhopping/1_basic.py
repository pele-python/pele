"""
Example 1: Simple basin hopping
"""
from pele.systems import LJCluster

natoms = 12
niter = 100
system = LJCluster(natoms)

db = system.create_database()
bh = system.get_basinhopping(database=db)
bh.run(niter)
print "the lowest energy found after", niter, " basinhopping steps is", db.minima()[0].energy
print ""

# some visualization
try: 
    import pele.utils.pymolwrapper as pym
    pym.start()
    pym.draw_spheres(bh.coords, "A", 1)
except:
    print "Could not draw using pymol, skipping this step"
