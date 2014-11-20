"""
Example 3: Saving all minima found to an xyz file
"""
from pele.systems import LJCluster
from pele.utils.xyz import write_xyz

natoms = 12
niter = 100
system = LJCluster(natoms)

db = system.create_database()
bh = system.get_basinhopping(database=db)
bh.run(niter)

with open("lowest", "w") as fout:
    for minimum in db.minima():
        title = "energy = ", str(minimum.energy)
        write_xyz(fout, minimum.coords, title)
           
############################################################
# some visualization
############################################################
try: 
    import pele.utils.pymolwrapper as pym
    pym.start()
    frame=1  
    for minimum in db.minima():        
        pym.draw_spheres(minimum.coords.reshape(-1, 3), "A", frame)
        frame += 1
except:
    print "Could not draw using pymol, skipping this step" 
