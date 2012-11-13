import numpy as np
from pygmin import takestep
from math import pi
from pygmin.utils.rbtools import CoordsAdapter
from pygmin.utils import rotations

# This is the takestep routine for OXDNA. It is a standard rigid body takestep
# routine, but I put it here to be able to start modifying it
class OXDNATakestep(takestep.TakestepInterface):
    def __init__(self, displace=1.0, rotate=0.5*pi):
        self.displace = displace
        self.rotate = rotate
        
    def takeStep(self, coords, **kwargs):
        # easy access to coordinates
        ca = CoordsAdapter(nrigid=coords.size/6, coords = coords)
        
        # random displacement for positions
        ca.posRigid[:] += 2.*self.displace*(np.random.random(ca.posRigid.shape)-0.5)
        
        # random rotation for angle-axis vectors
        takestep.rotate(self.rotate, ca.rotRigid)
        
    # this is necessary for adaptive step taking
    def scale(self, factor):
        self.rotate *= factor
        self.displace *= factor

    @property
    def stepsize(self):
        return [self.rotate, self.displace]
            
class OXDNAScrewStep(takestep.TakestepInterface):
    def __init__(self, rotate_backbone=0.5*pi, rotate_base=0., ntorsionmoves=1):
        self.rotate_backbone = rotate_backbone
        self.rotate_base = rotate_base
        self.ntorsionmoves=1

    def takeStep(self, coords, **kwargs):
        # easy access to coordinates
        ca = CoordsAdapter(nrigid=coords.size/6, coords = coords)

        for i in xrange(ca.nrigid):
            a = np.dot(rotations.aa2mx(ca.rotRigid[i]), np.array([1., 0., 0.]))
            a *= 2.*(np.random.random()-0.5)*self.rotate_base
            ca.rotRigid[i] = rotations.rotate_aa(ca.rotRigid[i], a)
               
         # random rotation for angle-axis vectors
        if self.rotate_backbone != 0.:
            for j in xrange(self.ntorsionmoves):
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
                p *= 2.*(np.random.random()-0.5)*self.rotate_backbone
                # convert random rotation to a matrix
                mx = rotations.aa2mx(p)
                # center of rotation is in middle of backbone bond
                center = 0.5*(x1 + x2)
                
                # apply rotation to positions and orientations
                for i in xrange(index, ca.nrigid):
                    a = np.dot(rotations.aa2mx(ca.rotRigid[i]), np.array([1., 0., 0.]))
                    ca.rotRigid[i] = rotations.rotate_aa(ca.rotRigid[i], p)
                    x = ca.posRigid[i] - 0.4*a
                    ca.posRigid[i] = np.dot(mx, x - center) + center
                    a = np.dot(rotations.aa2mx(ca.rotRigid[i]), np.array([1., 0., 0.]))
                    ca.posRigid[i]+=0.4*a
                
# this is necessary for adaptive step taking
    def scale(self, factor):
        self.rotate_backbone *= factor
        self.rotate_base *= factor

    @property
    def stepsize(self):
        return [self.rotate_backbone, self.rotate_base]
    
# this class should generate a fully random configuration
class OXDNAReseed(takestep.TakestepInterface):
    def __init__(self, radius=3.0):
        self.radius = radius
    
    def takeStep(self, coords, **kwargs):
        # easy access to coordinates
        ca = CoordsAdapter(nrigid=coords.size/6, coords = coords)
        
        # random displacement for positions
        ca.posRigid[:] = 2.*self.radius*(np.random.random(ca.posRigid.shape)-0.5)
        
        # random rotation for angle-axis vectors
        for rot in ca.rotRigid:
            rot[:] = rotations.random_aa()