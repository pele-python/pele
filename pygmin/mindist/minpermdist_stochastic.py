import numpy as np
from exact_match import StandardClusterAlignment
from pygmin.utils import rotations     
from _minpermdist_policies import TransformAtomicCluster, MeasureAtomicCluster

__all__ = ["MinPermDistCluster"]

class MinPermDistCluster(object):
    def __init__(self, niter=100, verbose=False, tol=0.01, accuracy=0.01,
                 measure=MeasureAtomicCluster(), transform=TransformAtomicCluster()):
        
        self.niter = 100
        
        self.verbose = verbose
        self.measure = measure
        self.transform=transform
        self.accuracy = accuracy
        self.tol = tol
        
    def check_match(self, x1, x2, rot, invert):
        x2_trial = x2.copy()
        if(invert):
            self.transform.invert(x2_trial)
        self.transform.rotate(x2_trial, rot)

        
        # get the best permutation
        dist, perm = self.measure.find_permutation(x1, x2_trial)
        x2_trial = self.transform.permute(x2_trial, perm)
       
        # now find best rotational alignment, this is more reliable than just
        # aligning the 2 reference atoms
        dist, rot2 = self.measure.find_rotation(x1, x2_trial)
        self.transform.rotate(x2_trial, rot2)
        # use the maximum distance, not rms as cutoff criterion
        
        dist =  self.measure.get_dist(x1, x2_trial)
        
        if dist < self.distbest:
            self.distbest = dist
            self.rotbest = np.dot(rot2, rot)
            self.x2_best = x2_trial    
    
    def finalize_best_match(self, x1):
        self.transform.translate(self.x2_best, self.com_shift)

        dist = self.measure.get_dist(x1, self.x2_best)
        if np.abs(dist - self.distbest) > 1e-6:
            raise RuntimeError        
        if self.verbose:
            print "finaldist", dist, "distmin", self.distbest

        return dist, self.x2_best
        
    def __call__(self, coords1, coords2):        
        # we don't want to change the given coordinates
        check_inversion = False
        x1 = np.copy(coords1)
        x2 = np.copy(coords2)
    
        com1 = self.measure.get_com(x1)
        self.transform.translate(x1, -com1)
        com2 = self.measure.get_com(x2)
        self.transform.translate(x2, -com2)

        self.com_shift = com1
        
        self.mxbest = np.identity(3)
        self.distbest = self.measure.get_dist(x1, x2)
    
        if self.distbest < self.tol:
            return self.distbest, x1, x2
        
#        for rot, invert in StandardClusterAlignment(x1, x2, accuracy=self.accuracy, 
#                                                    can_invert=self.transform.can_invert()):
        for rot, invert in StandardClusterAlignment(x1, x2):
            pass
            #self.check_match(x1, x2, rot, invert)
            #if self.distbest < self.tol:
            #    dist, x2 = self.finalize_best_match(coords1)
            #    return dist, coords1, x2
        
        # if we didn't find a perfect match here, try random rotations to optimize the match
        for i in range(self.niter):
            rot = rotations.aa2mx(rotations.random_aa())
            self.check_match(x1, x2, rot, False)
            if(self.transform.can_invert):
                self.check_match(x1, x2, rot, True)

#        self.transform.rotate(X2, mxbest)
#        dist, perm = self.measure.find_permutation(X1, X2)
#        X2 = self.transform.permute(X2, perm)
#        tmp, mx = self.measure.find_rotation(X1.copy(), X2.copy())
#        self.transform.rotate(X2, mx)
        
        # TODO: should we do an additional sanity check for permutation / rotation?        
        
        dist, x2 = self.finalize_best_match(coords1)                
        return dist, coords1, x2
    
def test(X1, X2, lj, atomtypes=["LA"], fname = "lj.xyz",
         minPermDist=MinPermDistCluster()):
    import copy
    natoms = len(X1) / 3
        
    X1i = copy.copy(X1)
    X2i = copy.copy(X2)
    
    printlist = []
    printlist.append((X2.copy(), "X2 initial"))
    printlist.append((X1.copy(), "X1 initial"))


    distinit = np.linalg.norm(X1-X2)
    print "distinit", distinit

    (dist, X1, X2) = minPermDist(X1,X2)
    distfinal = np.linalg.norm(X1-X2)
    print "dist returned    ", dist
    print "dist from coords ", distfinal
    print "initial energies (post quench)", lj.getEnergy(X1i), lj.getEnergy(X2i)
    print "final energies                ", lj.getEnergy(X1), lj.getEnergy(X2)

    printlist.append((X1.copy(), "X1 final"))
    printlist.append((X2.copy(), "X2 final"))


    import pygmin.printing.print_atoms_xyz as printxyz
    with open(fname, "w") as fout:
        for xyz, line2 in printlist:
            printxyz.printAtomsXYZ(fout, xyz, line2=line2 +" "+ str(lj.getEnergy(xyz)))
            
def test_LJ(natoms = 12, **kwargs):
    from pygmin.potentials.lj import LJ
    import pygmin.defaults
    import pygmin.utils.rotations as rot
    from pygmin.mindist.permutational_alignment import permuteArray
    import random
    
    quench = pygmin.defaults.quenchRoutine
    lj = LJ()
    X1 = np.random.uniform(-1,1,[natoms*3])*(float(natoms))**(1./3)
    #quench X1
    ret = quench( X1, lj.getEnergyGradient)
    X1 = ret[0]
    X2 = np.random.uniform(-1,1,[natoms*3])*(float(natoms))**(1./3)
    #make X2 a rotation of X1
    print "testing with", natoms, "atoms, with X2 a rotated and permuted isomer of X1"
    aa = rot.random_aa()
    rot_mx = rot.aa2mx( aa )
    for j in range(natoms):
        i = 3*j
        X2[i:i+3] = np.dot( rot_mx, X1[i:i+3] )
    perm = range(natoms)
    random.shuffle( perm )
    print perm
    X2 = permuteArray( X2, perm)

    #X1 = np.array( [ 0., 0., 0., 1., 0., 0., 0., 0., 1.,] )
    #X2 = np.array( [ 0., 0., 0., 1., 0., 0., 0., 1., 0.,] )
    import copy
    X1i = copy.copy(X1)
    X2i = copy.copy(X2)
 
    print "******************************"
    print "testing normal LJ  ISOMER"
    print "******************************"
    test(X1, X2, lj, **kwargs)
    
    print "******************************"
    print "testing normal LJ  non isomer"
    print "******************************"
    X2 = np.random.uniform(-1,1,[natoms*3])*(float(natoms))**(1./3)
    ret = quench( X2, lj.getEnergyGradient)
    X2 = ret[0]
    
    Y = X1.reshape([-1,3])
    Y+=np.random.random(3)
    X1[:] = Y.flatten()
 
    test(X1, X2, lj, **kwargs)
    

    distinit = np.linalg.norm(X1-X2)
    print "distinit", distinit

if __name__ == "__main__":
    test_LJ()