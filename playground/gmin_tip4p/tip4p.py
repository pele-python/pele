from math import sin, cos, pi
import numpy as np
from pele.angleaxis import RigidFragment, RBSystem
from pele.mindist.rmsfit import findrotation_kabsch
from pele.utils import rotations
from pele.transition_states import NEB, InterpolatedPath

def dump_path(filename, system, path):
    fl = open(filename, "w")
    lbls = system.get_atom_labels()
    for c in path:
        from pele.printing import printAtomsXYZ
        atomistic = system.to_atomistic(c)
        fl.write("%d\n\n"%len(lbls))
        for lbl, x in zip(lbls, atomistic):
            fl.write("%s %f %f %f\n"%(lbl, x[0], x[1], x[2]))
    fl.close()

def water():
    water = RigidFragment()
    rho   = 0.9572
    theta = 104.52/180.0*pi      
    water.add_atom("O", np.array([0., 0., 0.]), 16.)
    water.add_atom("H", rho*np.array([0.0, sin(0.5*theta), cos(0.5*theta)]), 1.)
    water.add_atom("H", rho*np.array([0.0, -sin(0.5*theta), cos(0.5*theta)]), 1.)
    water.finalize_setup()
    return water

def align(system, coords1, coords2):
    c1 = system.coords_adapter(coords1)
    c2 = system.coords_adapter(coords2)
    R = findrotation_kabsch(c2.posRigid, c1.posRigid).transpose()
    #R = rotations.aa2mx(p)
    for x, p in zip(c2.posRigid, c2.rotRigid):
        x[:] = np.dot(R, x)
        p[:] = rotations.rotate_aa(p, rotations.mx2aa(R))
    
    # now account for symmetry in water
    for p1, p2 in zip(c1.rotRigid,c2.rotRigid):
        theta1 = np.linalg.norm(rotations.rotate_aa(p2,-p1))
        p2n = rotations.rotate_aa(np.array([0., 0., pi]), p2)
        theta2 = np.linalg.norm(rotations.rotate_aa(p2n,-p1))
        theta1 -= int(theta1/2./pi)*2.*pi
        theta2 -= int(theta2/2./pi)*2.*pi
        if(theta2 < theta1): 
            p2[:]=p2n
            
def get_path(system, coords1, coords2, nimages):
    align(system, coords1, coords2)
    path=[x for x in InterpolatedPath(coords2, coords1, nimages, interpolator=system.interpolate)]
    system.align_path(path)
    return path

def path_dist_cart(system, path):
    cart = []
    for i in xrange(0,len(path)-1):
        cart.append(np.sqrt(np.linalg.norm(path[i] - path[i+1])**2))
    return cart
    
def path_dist_aa(system, path):
    aa = []
    for i in xrange(0,len(path)-1):
        aa.append(np.sqrt(system.distance_squared(path[i], path[i+1])))    
    return aa

