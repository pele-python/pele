import numpy as np
from exact_match import StandardClusterAlignment
from pele.utils import rotations     
from _minpermdist_policies import TransformAtomicCluster, MeasureAtomicCluster
from pele.mindist.periodic_exact_match import MeasurePeriodic,\
    TransformPeriodic, MeasurePeriodicRigid, TransformPeriodicRigid
from pele.utils.rbtools import CoordsAdapter
from time import ctime

__all__ = ["MinPermDistCluster"]

class MinPermDistCluster(object):
    """
    Minimize the distance between two clusters.  
    
    Parameters
    ----------
    niter : int
        the number of basinhopping iterations to perform
    verbose : boolean 
        whether to print status information
    accuracy :float, optional
        accuracy for standard alignment which determines if the structures are identical
    tol : float, optional
        tolerance for an exact match to stop iterations
    transform : 
        Transform policy which tells MinpermDist how to transform the given coordinates
    measure : 
        measure policy which tells minpermdist how to perform certains measures on the coordinates.
    
    Notes
    -----

    The following symmetries will be accounted for::
    
    1. Translational symmetry
    #. Global rotational symmetry
    #. Permutational symmetry
    #. Point inversion symmetry

    
    The algorithm here to find the best distance is
    
    for rotation in standardalignments:
        optimize permutation
        optimize rotation
        check_match
        
    for i in range(niter):    
        random_rotation
        optimize permutations
        align rotation
        check_match
        
    The minpermdist algorithm is generic and can act on various types of
    coordinates, e.g. carthesian, angle axis, .... The transform and measure
    policies define and interface to manipulate and analyze a given set of
    coordinates. If the coordinates don't have a standard format, custom policies
    can be specified. As an example see the angle axis minpermdist routines.
        
    See also
    --------
    TransformPolicy, MeasurePolicy
    
    """
    def __init__(self, niter=10, verbose=False, tol=0.01, accuracy=0.01,
                 measure=MeasureAtomicCluster(), transform=TransformAtomicCluster()):
        
        self.niter = niter
        
        self.verbose = verbose
        self.measure = measure
        self.transform=transform
        self.accuracy = accuracy
        self.tol = tol
        
    def check_match(self, x1, x2, rot, invert):
        ''' check a given rotation for a match '''
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
            self.invbest = invert
            self.x2_best = x2_trial    
    
    def finalize_best_match(self, x1):
        ''' do final processing of the best match '''
        self.transform.translate(self.x2_best, self.com_shift)

        dist = self.measure.get_dist(x1, self.x2_best)
        if np.abs(dist - self.distbest) > 1e-6:
            raise RuntimeError        
        if self.verbose:
            print "finaldist", dist, "distmin", self.distbest

        return dist, self.x2_best

    def _standard_alignments(self, x1, x2):
        ''' get iterator for standard alignments '''
        return StandardClusterAlignment(x1, x2, accuracy=self.accuracy, 
                                        can_invert=self.transform.can_invert())
        
    def __call__(self, coords1, coords2):        
        '''
        Parameters
        ----------
        coords1, coords2 : np.array 
            the structures to align.  X2 will be aligned with X1, both
            the center of masses will be shifted to the origin
            
        Returns
        -------
        a triple of (dist, coords1, coords2). coords1 are the unchanged coords1
        and coords2 are brought in best alignment with coords2
        '''
        # we don't want to change the given coordinates
        coords1 = coords1.copy()
        coords2 = coords2.copy()
        
        x1 = np.copy(coords1)
        x2 = np.copy(coords2)
    
        com1 = self.measure.get_com(x1)
        self.transform.translate(x1, -com1)
        com2 = self.measure.get_com(x2)
        self.transform.translate(x2, -com2)

        self.com_shift = com1
        
        self.mxbest = np.identity(3)
        self.distbest = self.measure.get_dist(x1, x2)
        self.x2_best = x2.copy()
        
        # sn402: The unlikely event that the structures are already nearly perfectly aligned.
        if self.distbest < self.tol:
            dist, x2 = self.finalize_best_match(coords1)
            return self.distbest, coords1, x2
        
        # sn402: This gives a systematic search for more-or-less perfect matches between x1 and x2.
        for rot, invert in self._standard_alignments(x1, x2):
            self.check_match(x1, x2, rot, invert)
            if self.distbest < self.tol:
                dist, x2 = self.finalize_best_match(coords1)
                return dist, coords1, x2
        
        # if we didn't find a perfect match here, try random rotations to optimize the match
        for i in range(self.niter):
            rot = rotations.aa2mx(rotations.random_aa())
            self.check_match(x1, x2, rot, False)
            if(self.transform.can_invert()):
                self.check_match(x1, x2, rot, True)

        # TODO: should we do an additional sanity check for permutation / rotation?        
        
        dist, x2 = self.finalize_best_match(coords1)                
        return dist, coords1, x2
    
    
class MinPermDistBulk(object):
    def __init__(self, boxvec, measure, transform=TransformPeriodicRigid(), niter=10, verbose=False, tol=0.01, 
                 accuracy=0.01):        
        self.niter = niter       
        self.verbose = verbose
        self.measure = measure
        self.transform=transform
        self.accuracy = accuracy
        self.tol = tol
        self.boxvec = boxvec
        
    def __call__(self, coords1, coords2):        
        '''
        Parameters
        ----------
        coords1, coords2 : np.array 
            the structures to align.  X2 will be aligned with X1
            Both structures are arrays of rigid body com positions and aa vectors
            
        Returns
        -------
        a triple of (dist, coords1, coords2). coords1 are the unchanged coords1
        and coords2 are brought in best alignment with coords2
        '''
        # we don't want to change the given coordinates
        coords1 = coords1.copy()
        coords2 = coords2.copy()
        
        x1 = np.copy(coords1)
        x2 = np.copy(coords2)
        self.distbest = self.measure.get_dist(x1, x2)
        
        ca1 = CoordsAdapter(coords=coords1)
        ca2 = CoordsAdapter(coords=coords2)        

        dx = ca2.posRigid - ca1.posRigid
        dx -= np.round(dx / self.boxvec) * self.boxvec
        ave = dx.sum(0)/ca1.nrigid
        
        self.transform.translate(x2, -ave)
        self.distbest = self.measure.get_dist(x1,x2)
     
        dist, x2 = self.finalize_best_match(coords1, x2)          
        return dist, coords1, x2        

    def finalize_best_match(self, x1, best_x2):
        ''' do final processing of the best match '''
        ca = CoordsAdapter(coords=best_x2)
        ca.posRigid -= np.round(ca.posRigid / self.boxvec) * self.boxvec

        dist = self.measure.get_dist(x1, best_x2)
        #if np.abs(dist - self.distbest) > 1e-6:
        #    raise RuntimeError        
        if self.verbose:
            print "finaldist", dist, "distmin", self.distbest
        return dist, best_x2

#
# testing only below here
#

def test(X1, X2, lj, atomtypes=["LA"], fname = "lj.xyz",
         minPermDist=MinPermDistCluster()): # pragma: no cover
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


    import pele.printing.print_atoms_xyz as printxyz
    with open(fname, "w") as fout:
        for xyz, line2 in printlist:
            printxyz.printAtomsXYZ(fout, xyz, line2=line2 +" "+ str(lj.getEnergy(xyz)))
            
def test_LJ(natoms = 12, **kwargs): # pragma: no cover
    from pele.potentials.lj import LJ
    from pele.optimize import mylbfgs
    import pele.utils.rotations as rot
    from pele.mindist.permutational_alignment import permuteArray
    import random
    
    quench = mylbfgs
    lj = LJ()
    X1 = np.random.uniform(-1,1,[natoms*3])*(float(natoms))**(1./3)
    # quench X1
    ret = quench( X1, lj)
    X1 = ret.coords
    X2 = np.random.uniform(-1,1,[natoms*3])*(float(natoms))**(1./3)
    # make X2 a rotation of X1
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
    ret = quench( X2, lj)
    X2 = ret.coords
    
    Y = X1.reshape([-1,3])
    Y+=np.random.random(3)
    X1[:] = Y.flatten()
 
    test(X1, X2, lj, **kwargs)
    

    distinit = np.linalg.norm(X1-X2)
    print "distinit", distinit


def sn402_test():
    import time
    start = time.time()
    np.random.seed(0)
    boxvec = np.array([10.,10.,10.])
    import pele.angleaxis._otp_bulk as OTP
    system = OTP.OTPBulk(2,boxvec,2.5)
    
    thismeasure = MeasurePeriodicRigid(box_lengths=boxvec, topology=system.aatopology)

    a = MinPermDistBulk(boxvec, thismeasure)

    #X1 = np.array([2.,1.,2.,2.,4.,5.,0.,0.,0.,0.,0.,0.])
    #X2 = np.array([5.,1.,2.,3.,4.,5.,0.,0.,0.,0.,0.,0.])
    X1 = np.random.random(36)*10-5
    #displacement = np.random.random(18)*5-5
    displacement = np.array([8,0,0,8,0,0,8,0,0,8,0,0,8,0,0,8,0,0,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    X2 = X1+displacement
    dist, x1, newx2 = a(X1, X2)
    
    end = time.time()
    print dist
    print x1
    print displacement
    print newx2-X2
    print newx2
    print "time elapsed:", end-start
    
    
if __name__ == "__main__":
    #test_LJ()
    sn402_test()