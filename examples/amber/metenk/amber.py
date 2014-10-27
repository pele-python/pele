import ambgmin_ as GMIN
import pele.potentials.gminpotential as gminpot
import numpy as np
import pele.basinhopping as bh
from pele.optimize import _quench as quench
from pele.takestep import displace

# export PYTHONPATH=/home/ss2029/svn/GMIN/bin:$PWD/../..

GMIN.initialize()   
pot = gminpot.GMINPotental(GMIN)

coords = pot.getCoords()

step = displace.RandomDisplacement(stepsize=0.7)

opt = bh.BasinHopping(coords, pot, takeStep=step, quench=quench.lbfgs_py)
opt.quenchParameters['tol'] = 1e-4
opt.run(3)

# some visualization
try: 
    import pele.utils.pymolwrapper as pym
    pym.start()
    pym.draw_spheres(opt.coords, "A", 1)
except:
    print "Could not draw using pymol, skipping this step"
from pele.utils.xyz import write_xyz
write_xyz(open("final.xyz", "w"), opt.coords)
