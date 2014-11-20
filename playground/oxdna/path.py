import numpy as np
from pele.utils.rbtools import CoordsAdapter
from pele.utils import rotations
import oxdnagmin_ as GMIN
from pele.potentials import GMINPotential
from pele.transition_states import NEB, InterpolatedPath
from pele.storage import Database
from pele.mindist.rmsfit import findrotation_kabsch
from math import pi
from pele.angleaxis import aamindist
from oxgui import OXDNASystem

infile="path/int.2min.xyz"

def map_to_aa(xyz):
    coords = np.zeros(6*13)
    ca = CoordsAdapter(nrigid=13, coords=coords)
    for i in xrange(13):
        ca.posRigid[i] = 0.5*(xyz[2*i] + xyz[2*i+1])    
        a1 = -(xyz[2*i] - xyz[2*i+1])
        
        if(i<12):
            a3 = -(xyz[2*i+1] - xyz[2*i+3])
        else:
            a3 = -(xyz[2*i-1] - xyz[2*i+1])
            
        a2 = -np.cross(a1, a3)
        a3 = np.cross(a1, a2)
        a1/=np.linalg.norm(a1)
        a2/=np.linalg.norm(a2)
        a3/=np.linalg.norm(a3)
 #       print np.dot(a1, a2), np.dot(a1, a3), np.dot(a2, a3)
        mx = np.array([a1, a2, a3]).transpose()
        p = rotations.mx2aa(mx)
        #print mx - rotations.aa2mx(p)
#        print a1-np.dot(rotations.aa2mx(p), np.array([1., 0., 0.]))
        ca.rotRigid[i] = p
    return coords
    
def export_xyz(fl, coords):
    ca = CoordsAdapter(nrigid=coords.size/6, coords = coords)
    fl.write("%d\n\n"%(2*ca.nrigid))
    for i in xrange(ca.nrigid):
        a = np.dot(rotations.aa2mx(ca.rotRigid[i]), np.array([1., 0., 0.]))
        x_back = ca.posRigid[i] - 0.4*a # backbone bead
        x_stack = ca.posRigid[i] + 0.4*a
        
        a = np.dot(rotations.aa2mx(ca.rotRigid[i]), np.array([0., 0., 1.]))
        x_tmp = x_back + 0.2*a
        
        fl.write("C %f %f %f\n"%(x_back[0], x_back[1], x_back[2]))
        fl.write("H %f %f %f\n"%(x_stack[0], x_stack[1], x_stack[2]))
        
system = OXDNASystem()
pot = system.get_potential()

#fin = open("test.xyz") #path/int.102min.xyz")
fin = open(infile)

path_xyz = []

while True:
    if(fin.readline() == ""):
        break
    
    fin.readline()
    coords = np.zeros([26,3])
    for i in xrange(26):
        x = [ float(x) for x in fin.readline().split()[1:]]
        coords[i] = x
    path_xyz.append(coords)
    
print "length of path:", len(path_xyz)

path = [ map_to_aa(xyz) for xyz in path_xyz]

# this block aligns the images
transform = aamindist.TransformAngleAxisCluster(system.aasystem)
measure = aamindist.MeasureAngleAxisCluster(system.aasystem, transform=transform)

com = measure.get_com(path[0])
transform.translate(path[0], -com)
for i in xrange(1,len(path)):
    com = measure.get_com(path[i])
    transform.translate(path[i], -com)
    dist, rot = measure.find_rotation(path[i-1], path[i])
    transform.rotate(path[i], rot)

#p#ath = [ x for x in IntterpolatedPath(db.minima[19], )]
#path = [ x for x in InterpolatedPath(path[0].copy(), path[-1].copy(), 34) ]
traj = open("traj.xyz", "w")
#for x in path:
#    #export_xyz(traj, x)
#    #ret = quench.mylbfgs(x, pot.getEnergyGradient)
#    #print i,pot.getEnergy(x), ret[1]
#    
#    # export_xyz(traj, ret[0])
for x in path:
    export_xyz(traj, x)

import pickle
pickle.dump(path, open("interpolate.pickle", "w"))
exit()
db=Database(db="oxdna.sqlite")
path[0]=db.minima()[19].coords
path[-1]=db.minima()[0].coords

e1 = []
e2 = []

e1.append(pot.getEnergy(path[0]))
e2.append(pot.getEnergy(path[0]))

#for i in xrange(1):
for i in xrange(len(path) - 1):
    e1.append(pot.getEnergy(path[i + 1]))
    c1 = CoordsAdapter(nrigid=13, coords=path[i])
    c2 = CoordsAdapter(nrigid=13, coords=path[i + 1])
    com1 = np.sum(c1.posRigid, axis=0) / float(13)
    com2 = np.sum(c1.posRigid, axis=0) / float(13)
    c1.posRigid -= com1
    c2.posRigid -= com2
    mx = findrotation_kabsch(c2.posRigid, c1.posRigid).transpose()
    # print mx
    c2.posRigid[:] = np.dot(mx, c2.posRigid.transpose()).transpose()
    for p in c2.rotRigid:
        p[:] = rotations.rotate_aa(p, rotations.mx2aa(mx))
    e2.append(pot.getEnergy(path[i + 1]))

    for p1, p2 in zip(c1.rotRigid, c2.rotRigid):
        n2 = p2 / np.linalg.norm(p2) * 2. * pi

        while True:
            p2n = p2 + n2
            if (np.linalg.norm(p2n - p1) > np.linalg.norm(p2 - p1)):
                break
            p2[:] = p2n

        while True:
            p2n = p2 - n2
            if (np.linalg.norm(p2n - p1) > np.linalg.norm(p2 - p1)):
                break
            p2[:]=p2n
            
    NEBquenchParams = dict()
    NEBquenchParams["nsteps"]=2000
    NEBquenchParams["maxErise"]=1e-1
    NEBquenchParams["tol"]=1e-5
    NEBquenchParams["iprint"]=1
    NEBquenchParams["maxstep"]=.1
    neb = NEB(path, pot, k=10., dneb=True, with_springenergy=True, quenchParams=NEBquenchParams)
    neb.optimize()

    path = neb.coords


print np.linalg.norm(path[1] - path[0])

for x in neb.coords:
    export_xyz(traj, x)
    
#np.savetxt("energies.txt", neb.energies)
import pylab as pl
pl.plot(neb.energies)
#pl.plot(e1)
#pl.plot(e2, "x")
pl.show()

traj.close()

    
    
