from __future__ import print_function
from __future__ import absolute_import
import numpy as np
from collections import namedtuple

from ._minpermdist_policies import TransformAtomicCluster, MeasureAtomicCluster
from . import rmsfit

__all__= ["StandardClusterAlignment", "ExactMatchCluster"]

class StandardClusterAlignment(object):
    """
    class to iterate over standard alignments for atomic clusters

    Quickly determines alignments of clusters which are possible exact matches.
    It uses atoms which are far away from the center to determine possible
    rotations. The algorithm does the following:

    1) Get 2 reference atoms from structure 1 which are farthest away from center
       and are not linear
    2) Determine candidates from structure 2 which are in same shell
       as reference atoms from structure 1 (+- accuracy)
    3) loop over all candidate combinations to determine
       orientation and check for match. Skip directly if angle of candidates
       does not match angle of reference atoms in structure 1.

    Parameters
    ----------
    coords1 : np.array
        first coordinates
    coords2 : np.array
        second coordinates
    accuracy : float
        accuracy of shell for atom candidates in standard alignment
    can_invert : boolean
        is an inversion possible?

    Examples
    --------

    >> for rot, invert in StandardClusterAlignment(X1, X2):
    >>     print "possible rotation:",rot,"inversion:",invert

    """
    def __init__(self, coords1, coords2, accuracy = 0.01, can_invert=True):
        x1 = coords1.reshape([-1,3]).copy()
        x2 = coords2.reshape([-1,3]).copy()

        self.accuracy = accuracy
        self.can_invert = can_invert

        # calculate distance of all atoms
        R1 = np.sqrt(np.sum(x1*x1, axis=1))
        R2 = np.sqrt(np.sum(x2*x2, axis=1))

        # at least 2 atoms are needed
        # get atom most outer atom

        # get 1. reference atom in configuration 1
        # use the atom with biggest distance to com
        idx_sorted = R1.argsort()
        idx1_1 = idx_sorted[-1]

        # find second atom which is not in a line
        cos_best = 99.00
        for idx1_2 in reversed(idx_sorted[0:-1]):
            # stop if angle is larger than threshold
            cos_theta1 = np.dot(x1[idx1_1], x1[idx1_2]) / \
                (np.linalg.norm(x1[idx1_1])*np.linalg.norm(x1[idx1_2]))

            # store the best match in case it is a almost linear molecule
            if np.abs(cos_theta1) < np.abs(cos_best):
                cos_best = cos_theta1
                idx1_2_best = idx1_2

            if np.abs(cos_theta1) < 0.9:
                break

        idx1_2 = idx1_2_best

        # do a very quick check if most distant atom from
        # center are within accuracy
        if np.abs(R1[idx1_1] - R2.max()) > accuracy:
            candidates1 = []
            candidates2 = []
        else:
            # get indices of atoms in shell of thickness 2*accuracy
            candidates1 = np.arange(len(R2))[ \
                 (R2 > R1[idx1_1] - accuracy)*(R2 < R1[idx1_1] + accuracy)]
            candidates2 = np.arange(len(R2))[ \
                 (R2 > R1[idx1_2] - accuracy)*(R2 < R1[idx1_2] + accuracy)]

        self.x1 = x1
        self.x2 = x2
        self.idx1_1 = idx1_1
        self.idx1_2 = idx1_2
        self.idx2_1 = None
        self.idx2_2 = None
        self.invert = False

        self.cos_theta1 = cos_theta1
        self.candidates2 = candidates2

        self.iter1 = iter(candidates1)
        self.iter2 = iter(self.candidates2)

    def __iter__(self):
        return self

    def __next__(self):
        # obtain first index for first call
        if self.idx2_1 is None:
            self.idx2_1 = next(self.iter1)

        # toggle inversion if inversion is possible
        if self.can_invert and self.invert == False and self.idx2_2 is not None:
            self.invert = True
        else:
            # determine next pair of indices
            self.invert = False
            # try to increment 2nd iterator
            try:
                self.idx2_2 = next(self.iter2)
            except StopIteration:
                # end of list, start over again
                self.iter2 = iter(self.candidates2)
                # and increment iter1
                self.idx2_1 = next(self.iter1)
                self.idx2_2 = None
                return next(self)

        if self.idx2_1 == self.idx2_2:
            return next(self)

        x1 = self.x1
        x2 = self.x2
        idx1_1 = self.idx1_1
        idx1_2 = self.idx1_2
        idx2_1 = self.idx2_1
        idx2_2 = self.idx2_2

        assert idx1_1 is not None
        assert idx1_2 is not None
        assert idx2_1 is not None
        assert idx2_2 is not None

        # we can immediately trash the match if angle does not match
        try:
            cos_theta2 = np.dot(x2[idx2_1], x2[idx2_2]) / \
                (np.linalg.norm(x2[idx2_1])*np.linalg.norm(x2[idx2_2]))
        except ValueError:
            raise
        if np.abs(cos_theta2 - self.cos_theta1) > 0.5:
            return next(self)

        mul = 1.0
        if self.invert:
            mul=-1.0

        # get rotation for current atom match candidates
        dist, rot = rmsfit.findrotation(
            x1[[idx1_1, idx1_2]], mul*x2[[idx2_1, idx2_2]], align_com=False)

        return rot, self.invert

    def next(self):
        return self.__next__()


class ClusterTransoformation(object):
    """an object that defines a transformation on a cluster"""
    translation = None
    rotation = None
    permutation = None
    invert = False
    


class ExactMatchCluster(object):
    """ Deterministic check if 2 clusters are a perfect match

        Determines quickly if 2 clusters are a perfect match. It uses
        check_standard_alignment_cluster to get possible orientations.



        Parameters
        ----------

        accuracy: float, optional
            maximum deviation of atoms to still consider cluster as a match

        check_inversion: boolean, optional
            check for inversion symmetry, default is True

        permlist: iteratable, optional
            list of allowed permutations. Default is None which means all
            particles can be permuted

        align_com: boolean, optional
            Flag if com should be removed before comparison

        Examples
        --------

        >>> x1 = np.random.random(3*natoms)
        >>> x2 = x1 + 1e-4*np.random.random(x1.shape)
        >>> matches = ExactClusterMatch(accuracy=1e-3)
        >>> if match(x1, x2):
        >>>     print "the two structures are identical

    """
    def __init__(self, tol = 0.01, accuracy=0.01, transform=TransformAtomicCluster(), measure=MeasureAtomicCluster()):
        self.accuracy = accuracy
        self.tol = tol
        self.transform = transform
        self.measure = measure

    def standard_alignments(self, coords1, coords2):
        """return an iterator over the standard alignments"""
        return StandardClusterAlignment(coords1, coords2, accuracy = self.accuracy,
                                   can_invert=self.transform.can_invert())


    def __call__(self, coords1, coords2, check_inversion=True):
        """return True if the structures are an exact mach, False otherwise"""
        alignment = self.find_transformation(coords1, coords2, check_inversion=check_inversion)
        return alignment is not None

    def find_transformation(self, coords1, coords2, check_inversion=True):
        """Return None if the two structures are different, else return the transformations that brings them into alignment

        Returns
        -------
        None:
            if the two structures are different
        namedtuple:
            with fields: rotation, permutation, inversion
        """
        self.exact_transformation = ClusterTransoformation()

        com1 = self.measure.get_com(coords1)
        x1 = coords1.copy()
        self.transform.translate(x1, -com1)

        com2 = self.measure.get_com(coords2)
        x2 = coords2.copy()
        self.transform.translate(x2, -com2)

        self.exact_transformation.translation = com1 - com2

        for rot, invert in self.standard_alignments(x1, x2):
            if invert and not check_inversion: continue
            ret = self.check_match(x1, x2, rot, invert)
            if ret is not None:
                self.exact_transformation.invert = invert
                self.exact_transformation.rotation = ret.rotation
                self.exact_transformation.permutation = ret.permutation
                return self.exact_transformation
        return None

    def check_match(self, x1, x2, rot, invert):
        """Give a rotation and inversion, make a more detailed comparison if the 2 structures match

        Parameters
        ----------

        rot: np.array, 3x3
            guessed rotation based on reference atoms
        invert: boolean
            True do match for inverted coordinates

        returns: boolean
            True or False for match

        """
        x2_trial = x2.copy()

        # apply the inversion
        # note that inversion happens before rotation.  this is important as they are not commutative (I think)
        if invert:
            self.transform.invert(x2_trial)

        # apply the rotation to x2_trial
        self.transform.rotate(x2_trial, rot)

        # get the best permutation
        dist, perm = self.measure.find_permutation(x1, x2_trial)
        x2_trial = self.transform.permute(x2_trial, perm)

        # now find best rotational alignment, this is more reliable than just
        # aligning the 2 reference atoms
        dist, rot2 = self.measure.find_rotation(x1, x2_trial)
        # further rotate x2_trial by rot2
        self.transform.rotate(x2_trial, rot2)
        # use the maximum distance, not rms as cutoff criterion

        self._last_checked_rotation = np.dot(rot2, rot)
        if  self.measure.get_dist(x1, x2_trial) < self.tol:
            # return the rotation and permutation as a named tupple
            ret = namedtuple("alignment", ["permutation", "rotation"])
            ret.permutation = perm
            ret.rotation = np.dot(rot2, rot)
            return ret
        return None

    def apply_transformation(self, x, tform):
        # apply permutation
        if tform.permutation is not None:
            self.transform.permute(x, tform.permutation)

        # subtract center of mass
        com = self.measure.get_com(x)
        self.transform.translate(x, -com)

        # invert
        if tform.invert and self.transform.can_invert():
            self.transform.invert(x)

        # rotate
        if tform.rotation is not None:
            self.transform.rotate(x, tform.rotation)

        # translate
        if tform.translation is not None:
            self.transform.translate(x, tform.translation)

        # add back in center of mass
        self.transform.translate(x, com)

#
# only testing below here
#

def test(): # pragma: no cover
    natoms = 35
    from pele.utils import rotations

    for i in range(100):
        xx1 = np.random.random(3*natoms)*5
        xx1 = xx1.reshape([-1,3])
        mx = rotations.q2mx(rotations.random_q())
        xx2 = -np.dot(mx, xx1.transpose()).transpose()
        xx2 +=2.*(np.random.random(xx2.shape)-0.5)*0.001
        #xx2 = xx1.copy()
        tmp = xx2[1].copy()
        xx2[1] = xx2[4]
        xx2[4] = tmp
        print(i, ExactMatchCluster()(xx1.flatten(), xx2.flatten()))


if __name__ == '__main__':
    test()

