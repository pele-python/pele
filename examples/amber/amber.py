import ambgmin_ as GMIN
import potentials.gminpotential as gminpot
import numpy as np
import basinhopping as bh
from optimize import quench
from takestep import displace

GMIN.initialize()   
pot = gminpot.GMINPotental(GMIN)

coords = pot.getCoords()

step = displace.RandomDisplacement(stepsize=0.7)

opt = bh.BasinHopping(coords, pot, takeStep=step, quenchRoutine=quench.lbfgs_py)
opt.quenchParameters['tol'] = 1e-4
opt.run(3)

# some visualization
try: 
    import utils.pymolwrapper as pym
    pym.start()
    pym.draw_spheres(opt.coords, "A", 1)
except:
    print "Could not draw using pymol, skipping this step"
import printing.print_atoms_xyz as pr
pr.printAtomsXYZ(open("final.xyz", "w"), opt.coords)
