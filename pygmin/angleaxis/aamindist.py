import numpy as np
from pygmin.utils import rotations
from pygmin.mindist import ExactMatchCluster, MinPermDistCluster, StandardClusterAlignment
from pygmin.mindist import TransformPolicy, MeasurePolicy
from pygmin.mindist import findrotation, find_best_permutation
from math import pi

class TransformAngleAxisCluster(TransformPolicy):
    ''' transformation rules for atomic clusters '''
    def __init__(self, topology):
        self.topology = topology
        self._can_invert = True
        for s in topology.sites:
            if s.inversion is None:
                self._can_invert = False
                
    def translate(self, X, d):
        ca = self.topology.coords_adapter(X)
        if(ca.nrigid > 0):
            ca.posRigid += d

        if(ca.natoms > 0):
            ca.posAtom += d
        
    def rotate(self, X, mx):
        ca = self.topology.coords_adapter(X)
        if(ca.nrigid > 0):
            ca.posRigid[:] = np.dot(mx, ca.posRigid.transpose()).transpose()
            dp = rotations.mx2aa(mx)
            for p in ca.rotRigid:
                p[:] = rotations.rotate_aa(p, dp)

        if(ca.natoms > 0):
            ca.posAtom[:] = np.dot(mx, ca.posAtom.transpose()).transpose()
        
    def can_invert(self):
        return self._can_invert
    
    def invert(self, X):
        ca = self.topology.coords_adapter(X)
        ca.posRigid[:] = - ca.posRigid 
        for p, site in zip(ca.rotRigid, self.topology.sites):
            p[:] = rotations.rotate_aa(rotations.mx2aa(site.inversion), p)
    
    def permute(self, X, perm):
        Xnew = X.copy()
        ca = self.topology.coords_adapter(X)
        ca_new = self.topology.coords_adapter(Xnew)
        
        ca_new.posRigid[:] = ca.posRigid[perm]
        ca_new.rotRigid[:] = ca.rotRigid[perm]
        
        return Xnew
    
class MeasureAngleAxisCluster(MeasurePolicy):
    ''' measure rules for atomic clusters '''
    
    def __init__(self, topology, transform=None):
        self.topology = topology
        if transform == None:
            transform= TransformAngleAxisCluster(topology)
        self.transform = transform
        
    def get_com(self, X):
        ca = self.topology.coords_adapter(X)
        
        com = np.zeros(3)
        
        if ca.natoms > 0:
            raise NotImplementedError
        
        if ca.nrigid > 0:
            com = ca.posRigid.sum(0) / ca.nrigid
        
        return com

    def align(self, coords1, coords2):
        c1 = self.topology.coords_adapter(coords1)
        c2 = self.topology.coords_adapter(coords2)
        
        # now account for symmetry in water
        for p1, p2, site in zip(c1.rotRigid,c2.rotRigid, self.topology.sites):
            theta_min = 10.
            mx2 = rotations.aa2mx(p2)
            mx1 = rotations.aa2mx(p1).transpose()
            mx =  np.dot(mx1, mx2)
            for rot in site.symmetries:
                mx_diff = np.dot(mx, rot)
                theta = np.linalg.norm(rotations.mx2aa(mx_diff))
                                       
                theta -= int(theta/2./pi)*2.*pi
                if(theta < theta_min): 
                    theta_min = theta
                    rot_best = rot
            p2[:] = rotations.rotate_aa(rotations.mx2aa(rot_best), p2)
                    

    def get_dist(self, X1, X2):
        x1 = X1.copy()
        x2 = X2.copy()
        self.align(x1, x2)
        return self.topology.distance_squared(x1, x2)
    
    def find_permutation(self, X1, X2):
        ca1 = self.topology.coords_adapter(X1)
        ca2 = self.topology.coords_adapter(X2)
        
        return find_best_permutation(ca1.posRigid, ca2.posRigid)
    
    def find_rotation(self, X1, X2):
        ca1 = self.topology.coords_adapter(X1)        
        ca2 = self.topology.coords_adapter(X2)        
        if ca1.natoms > 0:
            raise NotImplementedError
        
        
        dist, mx = findrotation(ca1.posRigid.flatten(), ca2.posRigid.flatten())
        X2trans = X2.copy()
        self.transform.rotate(X2trans, mx)
        
        return self.get_dist(X1, X2trans), mx
    
class MeasureRigidBodyCluster(MeasureAngleAxisCluster):
    def get_dist(self, X1, X2):
        x1 = X1.copy()
        x2 = X2.copy()
        self.align(x1, x2)
        atom1 = self.topology.to_atomistic(x1)
        atom2 = self.topology.to_atomistic(x2)
        return np.linalg.norm(atom1-atom2)
    
class ExactMatchAACluster(ExactMatchCluster):
    def __init__(self, topology, transform=None, measure=None, **kwargs):
        self.topology = topology
        
        if transform is None:
            transform = TransformAngleAxisCluster(topology)
        
        if measure is None:
            measure = MeasureAngleAxisCluster(topology, transform=transform)
        
        ExactMatchCluster.__init__(self, transform=transform, measure=measure, **kwargs)
        
    def __call__(self, coords1, coords2):
        x1 = coords1.copy()
        x2 = coords2.copy()
        
        ca1 = self.topology.coords_adapter(x1)
        ca2 = self.topology.coords_adapter(x2)
        
        com1 = self.measure.get_com(coords1)
        self.transform.translate(x1, -com1)
        
        com2 = self.measure.get_com(coords2)
        self.transform.translate(x2, -com2)
        
        for rot, invert in StandardClusterAlignment(ca1.posRigid, ca2.posRigid, accuracy = self.accuracy,
                                   can_invert=self.transform.can_invert()):
            if self.check_match(x1, x2, rot, invert):
                return True
        return False
    
class MinPermDistAACluster(MinPermDistCluster):
    def __init__(self, topology, transform=None, measure=None, **kwargs):
        self.topology = topology
        
        if transform is None:
            transform = TransformAngleAxisCluster(topology)
        
        if measure is None:
            measure = MeasureAngleAxisCluster(topology, transform=transform)
        
        MinPermDistCluster.__init__(self, transform=transform, measure=measure, **kwargs)

    def _standard_alignments(self, x1, x2):
        ca1 = self.topology.coords_adapter(x1)
        ca2 = self.topology.coords_adapter(x2)
        return StandardClusterAlignment(ca1.posRigid, ca2.posRigid, accuracy=self.accuracy, 
                                        can_invert=self.transform.can_invert())

    def finalize_best_match(self, x1):
        self.transform.translate(self.x2_best, self.com_shift)
        self.measure.align(x1, self.x2_best)
        
        dist = self.measure.get_dist(x1, self.x2_best)
        if np.abs(dist - self.distbest) > 1e-6:
            raise RuntimeError        
        if self.verbose:
            print "finaldist", dist, "distmin", self.distbest

        return dist, self.x2_best