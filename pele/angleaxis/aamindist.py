from __future__ import print_function
from math import sqrt
from math import pi


import numpy as np

from pele.utils import rotations
from pele.mindist import ExactMatchCluster, MinPermDistCluster, StandardClusterAlignment
from pele.mindist import TransformPolicy, MeasurePolicy
from pele.mindist import findrotation, find_best_permutation
from pele.angleaxis import _cpp_aa


class TransformAngleAxisCluster(TransformPolicy):
    """transformation rules for angle axis clusters """
    def __init__(self, topology):
        self.topology = topology
        self._can_invert = True
        for s in topology.sites:
            if s.inversion is None:
                self._can_invert = False
        
        # js850> note: this makes the c++ topology mandatory.  If we want we can make it optional
        # by catching any failure and setting cpp_transform to None
        try:
            self.cpp_transform = _cpp_aa.cdefTransformAACluster(self.topology)
        except TypeError:
            self.cpp_transform = None

        
    def translate(self, X, d):
        """apply a translation"""
        ca = self.topology.coords_adapter(X)
        if ca.nrigid > 0:
            ca.posRigid += d

        if ca.natoms is not None and ca.natoms > 0:
            ca.posAtom += d
    
    def _rotate_python(self, X, mx):
        """the slow pythonic version of rotate"""
        ca = self.topology.coords_adapter(X)
        if ca.nrigid > 0:
            # rotate the center of mass positions
            ca.posRigid[:] = np.dot(ca.posRigid, mx.transpose())
            
            # rotate the angle axis rotations
            dp = rotations.mx2aa(mx)
            for p in ca.rotRigid:
                p[:] = rotations.rotate_aa(p, dp)

        if ca.natoms is not None and ca.natoms > 0:
            ca.posAtom[:] = np.dot(mx, ca.posAtom.transpose()).transpose()
    
    def rotate(self, X, mx):
        """rotate the com + angle-axis position X by the rotation matrix mx
        """
        if self.cpp_transform is not None:
            return self.cpp_transform.rotate(X, mx)
        else:
            return self._rotate_python(X, mx)
        
    def can_invert(self):
        return self._can_invert
    
    def invert(self, X):
        """invert the structure"""
        ca = self.topology.coords_adapter(X)
        ca.posRigid[:] = - ca.posRigid 
        for p, site in zip(ca.rotRigid, self.topology.sites):
            p[:] = rotations.rotate_aa(rotations.mx2aa(site.inversion), p)
    
    def permute(self, X, perm):
        """apply a permutation
        
        Paramters
        ---------
        X : np array
            The configuration
        perm : list of integers
            a list of integers giving the new order of the molecules. The length
        of perm must be the number of rigid bodies.
        """
        Xnew = X.copy()
        ca = self.topology.coords_adapter(X)
        ca_new = self.topology.coords_adapter(Xnew)
        
        # The following lines just re-order posRigid and rotRigid accordingly.
        ca_new.posRigid[:] = ca.posRigid[perm]
        ca_new.rotRigid[:] = ca.rotRigid[perm]
        
        return Xnew

    
class MeasureAngleAxisCluster(MeasurePolicy):
    """measure rules for angle axis clusters """
    
    def __init__(self, topology, transform=None, permlist=None):
        self.topology = topology
        if transform is None:
            transform= TransformAngleAxisCluster(topology)
        self.transform = transform
        self.permlist = permlist
        
        # js850> note: this makes the c++ topology mandatory.  If we want we can make it optional
        # by catching any failure and setting cpp_transform to None
        self.cpp_measure = _cpp_aa.cdefMeasureAngleAxisCluster(self.topology)
        
    def get_com(self, X):
        """return the center of mass"""
        ca = self.topology.coords_adapter(X)
        
        com = np.zeros(3)
        
        if ca.natoms is not None and ca.natoms > 0:
            raise NotImplementedError
        
        if ca.nrigid > 0:
            com = ca.posRigid.sum(0) / ca.nrigid
        # note: js850> This is treating all rigid bodies as if the have the same mass.  This
        # is probably a bug and should be updated.  However we might actually want the
        # center of geometry, so maybe we should add a new function get_cog().
        return com

    def _align_pythonic(self, coords1, coords2):
        """the slow pythonic version of align"""
        # note: this minimizes the angle-distance, but it might be better to 
        # minimize the atomistic distance.  These are not always the same
        c1 = self.topology.coords_adapter(coords1)
        c2 = self.topology.coords_adapter(coords2)
        
        # now account for inner-molecular symmetry
        for p1, p2, site in zip(c1.rotRigid,c2.rotRigid, self.topology.sites):
            theta_min = 10.
            mx2 = rotations.aa2mx(p2)
            mx1 = rotations.aa2mx(p1).transpose()
            mx =  np.dot(mx1, mx2)
            # find the symmetry operation which puts p2 into best alignment with p1
            for rot in site.symmetries:
                mx_diff = np.dot(mx, rot)
                # theta is the rotation angle between p1 and p2 after 
                # applying the symmetry operation rot to site2
                theta = np.linalg.norm(rotations.mx2aa(mx_diff))
                
                # remove any extra factors of 2*pi
                theta -= int(theta / (2.*pi)) * 2.*pi
                if theta < theta_min:
                    theta_min = theta
                    rot_best = rot
            p2[:] = rotations.rotate_aa(rotations.mx2aa(rot_best), p2)

    def align(self, coords1, coords2):
        """align the rotations so that the atomistic coordinates will be in best alignment"""
        if self.cpp_measure is not None:
            return self.cpp_measure.align(coords1, coords2)
        else:
            return self._align_pythonic(coords1, coords2)
                    

    def get_dist(self, x1, x2):
        """compute the distance between two configurations"""
        x1 = x1.copy()
        x2 = x2.copy()
        self.align(x1, x2)
        return sqrt(self.topology.distance_squared(x1, x2))
    
    def find_permutation(self, X1, X2):
        """find the rotation which minimizes the distance between the structures"""
        ca1 = self.topology.coords_adapter(X1)
        ca2 = self.topology.coords_adapter(X2)
        
        return find_best_permutation(ca1.posRigid, ca2.posRigid, permlist=self.permlist)
    
    def find_rotation(self, X1, X2):
        """find the rotation which minimizes the distance between the structures"""
        ca1 = self.topology.coords_adapter(X1)        
        ca2 = self.topology.coords_adapter(X2)        
        if ca1.natoms is not None and ca1.natoms > 0:
            raise NotImplementedError
        
        # align the center of mass coordinates
        dist, mx = findrotation(ca1.posRigid.ravel(), ca2.posRigid.ravel())
        X2trans = X2.copy()
        self.transform.rotate(X2trans, mx)
        
        # compute and return the distance between the rotated coordinates
        return self.get_dist(X1, X2trans), mx
     
            
class MeasureRigidBodyCluster(MeasureAngleAxisCluster):
    """perform measurements on clusters of rigid bodies"""
   
# js850> this is commented because it is unnecessary and assumes a non-periodic system
#        It's faster than the alternate implementation in python, but probably not in c++
#        The results should be identical
#    def get_dist(self, x1, x2):
#        """return the distance between two configurations"""
#        x1 = x1.copy()
#        x2 = x2.copy()
#        self.align(x1, x2)
#        atom1 = self.topology.to_atomistic(x1)
#        atom2 = self.topology.to_atomistic(x2)
#        return np.linalg.norm(atom1-atom2)

    
class ExactMatchAACluster(ExactMatchCluster):
    """test whether two structure are exactly the same"""
    def __init__(self, topology, transform=None, measure=None, **kwargs):
        self.topology = topology
        
        if transform is None:
            transform = TransformAngleAxisCluster(topology)
        
        if measure is None:
            measure = MeasureAngleAxisCluster(topology, transform=transform)
        
        ExactMatchCluster.__init__(self, transform=transform, measure=measure, **kwargs)

    def standard_alignments(self, coords1, coords2):
        ca1 = self.topology.coords_adapter(coords1)
        ca2 = self.topology.coords_adapter(coords2)
        return StandardClusterAlignment(ca1.posRigid, ca2.posRigid, accuracy=self.accuracy,
                                        can_invert=self.transform.can_invert())
            
class MinPermDistAACluster(MinPermDistCluster):
    """minimize the distance between two structures"""
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
            print("final dist", dist, "minimum dist", self.distbest)

        return dist, self.x2_best
        

