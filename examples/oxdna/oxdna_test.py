import numpy as np
import oxdnagmin_ as GMIN
from pygmin.potentials.gminpotential import GMINPotential
import pygmin.basinhopping as bh
from pygmin.takestep import displace
from pygmin.storage.database import Database
from pygmin import takestep
from pygmin.utils import rotations
from pygmin.utils.rbtools import CoordsAdapter
from math import pi
import os
from pygmin.printing import printAtomsXYZ
from pygmin.optimize.quench import lbfgs_py

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

class OXDNAScrewTakestep(takestep.TakestepInterface):
    def __init__(self, stepsize=0.5*pi):
        self.stepsize = stepsize

    def takeStep(self, coords, **kwargs):
        # easy access to coordinates
        ca = CoordsAdapter(nrigid=coords.size/6, coords = coords)

         # random rotation for angle-axis vectors
        for j in xrange(ca.nrigid):
            # choose bond to rotate around, index is first bead that changes
            index = np.random.randint(1, ca.nrigid)
            
            # determine backbone beads
            a1 = np.dot(rotations.aa2mx(ca.rotRigid[index-1]), np.array([1., 0., 0.]))
            a2 = np.dot(rotations.aa2mx(ca.rotRigid[index]), np.array([1., 0., 0.]))
            x1 = ca.posRigid[index-1] - 0.4*a1 # backbone bead
            x2 = ca.posRigid[index]   - 0.4*a2 # backbone bead
            
            # get bond vector as axis of rotation + random magnitude
            p = x2 - x1
            p /= np.linalg.norm(p)
            p *= np.random.random()*self.stepsize
            # convert random rotation to a matrix
            mx = rotations.aa2mx(p)
            # center of rotation is in middle of backbone bond
            center = 0.5*(x1 + x2)
            
            # apply rotation to positions and orientations
            for i in xrange(index, ca.nrigid):
                a = np.dot(rotations.aa2mx(ca.rotRigid[i]), np.array([1., 0., 0.]))
                ca.rotRigid[i] = rotations.rotate_aa(p, ca.rotRigid[i])
                x = ca.posRigid[i] - 0.4*a
                ca.posRigid[i] = np.dot(mx, x - center) + center
                a = np.dot(rotations.aa2mx(ca.rotRigid[i]), np.array([1., 0., 0.]))
                ca.posRigid[i]+=0.4*a
            
# this is necessary for adaptive step taking
    def scale(self, factor):
        self.stepsize *= factor
        
# initialize GMIN
GMIN.initialize()
# create a potential which calls GMIN
potential = GMINPotential(GMIN)
# get the initial coorinates
coords=potential.getCoords()
# create takestep routine

step = OXDNAScrewTakestep()
#os.system(TO_PDB%"before.dat")
out = open("traj.xyz", "w")
export_xyz(out, coords)
for i in xrange(1,10):
    step.takeStep(coords)
    coords_opt, E, rms, fcalls = lbfgs_py(coords, potential.getEnergyGradient)
    print E, fcalls
    export_xyz(out, coords)
    coords = coords_opt
    export_xyz(out, coords)
    

out.close()
