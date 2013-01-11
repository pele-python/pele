import numpy as np
from pygmin.mindist import ExactMatchCluster
from mindistutils import getAlignRotation
from permutational_alignment import findBestPermutation
from pygmin.utils import rotations

class NotImplemented(BaseException):
    pass

class TransformPolicy(object):   
    def translate(self, X, d):
        raise NotImplemented
    
    def rotate(self, X, mx):
        raise NotImplemented
    
    def can_invert(self):
        raise NotImplemented
    
    def invert(self, X):
        raise NotImplemented
    
class MeasurePolicy(object):
    def get_com(self, X):
        raise NotImplemented
    
    def get_dist(self, X1, X2):
        raise NotImplemented
    
    def find_permutation(self, X1, X2):
        ''' find the best permutation for the current set of coordinates '''
        raise NotImplemented
    
    def find_rotation(self, X1, X2):
        ''' find optimal rotation for 2 structures '''
        raise NotImplemented

class TransformAtomicCluster(TransformPolicy):
    # TODO: make this a class method?    
    def translate(self, X, d):
        ''' translate the system by a constant vector '''
        Xtmp = X.reshape([-1,3])
        Xtmp+=d
        
    def rotate(self, X, mx,):
        ''' rotate the system '''
        Xtmp = X.reshape([-1,3])
        Xtmp = np.dot(mx, Xtmp.transpose()).transpose()
        X[:] = Xtmp.reshape(X.shape)
        
    def can_invert(self):
        return True
    
    def invert(self, X):
        X[:] = -X
        
class MeasureAtomicCluster(MeasurePolicy):
    def get_com(self, X):
        ''' calculate the com '''
        X = np.reshape(X, [-1,3])
        natoms = len(X[:,0])
        com = X.sum(0) / natoms
        return com

    def get_dist(self, X1, X2):
        ''' calculate the distance between 2 sets of coordinates '''
        return np.linalg.norm(X1-X2)
    
    def find_permutation(self, X1, X2, permlist=None):
        ''' find the best permutation for the current set of coordinates '''
        return findBestPermutation(X1, X2, permlist)
    
    def find_rotation(self, X1, X2):
        ''' find optimal rotation for 2 structures '''
        dist, Q2 = getAlignRotation(X1, X2)
        mx = rotations.q2mx(Q2)
        return dist, mx     

class MinPermDistCluster(object):
    def __init__(self, niter=100,
                 measure=MeasureAtomicCluster(), transform=TransformAtomicCluster()):
        
        self.niter = 100
        
        self.verbose = False
        self.permlist = None
        self.measure = measure
        self.transform=transform
    
    def exact_match(self, X1, X2, accuracy = 0.01):
        ''' checks for an exact match '''
        exactmatch = ExactMatchCluster(accuracy=accuracy, permlist=self.permlist)
        return exactmatch(X1, X2)
    
    def _optimize_perm_rot(self, X1, X2):
        # TODO: add back option to quench again
        distbest = self.measure.get_dist(X1, X2)
        mxbest = np.identity(3)
        
        for i in range(self.niter):
            X2_ = X2.copy()
        
            #get and apply a random rotation
            aa = rotations.random_aa()
            mx = rotations.aa2mx(aa)
            mxtot = mx
            #print "X2.shape", X2.shape
            
            self.transform.rotate(X2_, mx)
            
            #optimize the permutations
            dist, X1, X2_ = self.measure.find_permutation(X1, X2_)
            if self.verbose:
                print "dist", dist, "distbest", distbest
            #print "X2.shape", X2.shape
            
            #optimize the rotation
            dist, mx2 = self.measure.find_rotation(X1, X2_)
            mxtot = np.dot(mx2, mxtot)
            
            # TODO: check with custom distance function here?
            # dist = self.get_dist(X1, X2)
            
            if dist < distbest:
                distbest = dist
                mxbest = mxtot
        return distbest, mxbest 
    
    def __call__(self, coords1, coords2, accuracy=0.01):
        # we don't want to change the given coordinates
        check_inversion = False
        X1 = np.copy(coords1)
        X2 = np.copy(coords2)
    
        #first check for exact match        
        if self.exact_match(X1, X2):
            #this is kind of cheating, I would prefer to return
            #X2 in best alignment and the actual (small) distance
            return 0.0, X1, X1.copy()
        
        com1 = self.measure.get_com(X1)
        self.transform.translate(X1, -com1)
        com2 = self.measure.get_com(X2)
        self.transform.translate(X2, -com2)
        
        #find the best rotation stochastically
        distbest, mxbest = self._optimize_perm_rot(X1, X2)
        
        use_inversion = False
        if self.transform.can_invert():
            X2i = X2.copy()
            self.transform.invert(X2i)
            distbest1, mxbest1 = self._optimize_perm_rot(X1, X2i)
            if distbest1 < distbest:
                if self.verbose:
                    print "using inversion in minpermdist"
                print "inversion rules", distbest1
                use_inversion = True
                distbest = distbest1
                mxbest = mxbest1
    
        #now we know the best rotation
        if use_inversion: X2 = X2i
        
        self.transform.rotate(X2, mxbest)
        dist, X1, X2 = self.measure.find_permutation(X1, X2)
        
        tmp, mx = self.measure.find_rotation(X1.copy(), X2.copy())
        self.transform.rotate(X2, mx)
        
        # Now perform a sanity check
        dist = self.measure.get_dist(X1, X2)
        if dist > distbest+0.001:
            print "ERROR: minPermDistRanRot: dist is different from distbest %f %f" % (dist, distbest)
        if self.verbose:
            print "finaldist", dist, "distmin", distbest
        
        # shift back com
        self.transform.translate(X1, com1)
        self.transform.translate(X2, com1)
        
        return dist, X1, X2
    
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