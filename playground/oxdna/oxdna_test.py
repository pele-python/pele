import numpy as np
import oxdnagmin_ as GMIN
from pele.potentials.gminpotential import GMINPotential
import pele.basinhopping as bh
from pele.takestep import displace
from pele.storage.database import Database
from pele import takestep
from pele.utils import rotations
from pele.utils.rbtools import CoordsAdapter
from math import pi
import os
from pele.printing import printAtomsXYZ
from pele.optimize.quench import lbfgs_py
from pele.systems.oxdna import OXDNATakestep

TO_PDB="python /home/vr274/opt/oxDNA/UTILS/traj2vis.py  pdb %s gmindnatop"

def export_xyz(fl, coords):
    ca = CoordsAdapter(nrigid=coords.size/6, coords = coords)
    fl.write("%d\n\n"%(2*ca.nrigid))
    for i in xrange(ca.nrigid):
        a = np.dot(rotations.aa2mx(ca.rotRigid[i]), np.array([1., 0., 0.]))
        x_back = ca.posRigid[i] - 0.4*a # backbone bead
        x_stack = ca.posRigid[i] + 0.4*a
        
        fl.write("C %f %f %f\n"%(x_back[0], x_back[1], x_back[2]))
        fl.write("H %f %f %f\n"%(x_stack[0], x_stack[1], x_stack[2]))
        
# initialize GMIN
GMIN.initialize()
# create a potential which calls GMIN
potential = GMINPotential(GMIN)
# get the initial coorinates
coords=potential.getCoords()
# create takestep routine

step = OXDNATakestep(displace=0, rotate_around_backbone=True)
#os.system(TO_PDB%"before.dat")
out = open("traj.xyz", "w")
export_xyz(out, coords)
for i in xrange(1,10):
    step.takeStep(coords)
    #coords_opt, E, rms, fcalls = lbfgs_py(coords, potential.getEnergyGradient)
    #print E, fcalls
    export_xyz(out, coords)
    #coords = coords_opt
    #export_xyz(out, coords)
    

out.close()
