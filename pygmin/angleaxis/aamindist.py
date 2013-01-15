import numpy as np
from pygmin.utils import rotations
from pygmin.mindist import TransformPolicy, MeasurePolicy
from pygmin.mindist import getAlignRotation, find_best_permutation

class TransformAngleAxisCluster(TransformPolicy):
    ''' transformation rules for atomic clusters '''
    def __init__(self, topology):
        self.topology = topology
        
    def translate(self, X, d):
        ca = self.topology.coords_adapter(X)
        if(ca.nrigid > 0):
            ca.posRigid += d

        if(ca.natoms > 0):
            ca.posAtom += d
        
    def rotate(self, X, mx,):
        ca = self.topology.coords_adapter(X)
        if(ca.nrigid > 0):
            ca.posRigid[:] = np.dot(mx, ca.posRigid.transpose()).transpose()
            dp = rotations.mx2aa(mx)
            for p in ca.rotRigid:
                p[:] = rotations.rotate_aa(p, dp)

        if(ca.natoms > 0):
            ca.posAtom[:] = np.dot(mx, ca.posAtom.transpose()).transpose()
        
    def can_invert(self):
        return False
    
    def invert(self, X):
        raise RuntimeError("system cannot be inverted")
    
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

    def get_dist(self, X1, X2):
        return self.topology.distance_squared(X1, X2)
    
    def find_permutation(self, X1, X2):
        ca1 = self.topology.coords_adapter(X1)
        ca2 = self.topology.coords_adapter(X2)
        
        return find_best_permutation(ca1.posRigid, ca2.posRigid)
    
    def find_rotation(self, X1, X2):
        ca1 = self.topology.coords_adapter(X1)        
        ca2 = self.topology.coords_adapter(X2)        
        if ca1.natoms > 0:
            raise NotImplementedError
        
        
        dist, Q2 = getAlignRotation(ca1.posRigid.flatten(), ca2.posRigid.flatten())
        mx = rotations.q2mx(Q2)
        X2trans = X2.copy()
        self.transform.rotate(X2trans, mx)
        
        return self.get_dist(X1, X2trans), mx
    