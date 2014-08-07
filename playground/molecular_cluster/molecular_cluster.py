'''
    permlists are messy; Jake hates them. but while we want for a better solution...
    for molecular clusters:
        if an element of permlist is a list of ints:
            permute all atoms within the list
        if an element of permlist is a list of lists:
            permute all atoms within the first list, then apply the best permutation
            to the atoms in subsequent lists
    e.g. three water molecules,
        permlist = [ [ [ O1 O2 O3 ], [ 1H1, 1H2, 1H3 ], [ 2H1, 2H2, 2H3 ] ],
                       [ 1H1, 2H1 ],
                       [ 1H2, 2H2 ],
                       [ 1H3, 2H3 ], ]
        where O1, etc. are integer indices, will permute all atoms associated with a
        pair of water molecules, then each pair of hydrogen atoms
    the routines here would work perfectly well for atomic clusters.
 '''

from pele.mindist._minpermdist_policies import MeasureAtomicCluster, TransformAtomicCluster
from pele.mindist.permutational_alignment import _make_cost_matrix, _find_permutations, find_permutations_hungarian
from pele.mindist import MinPermDistCluster, ExactMatchCluster
from pele.systems import AtomicCluster
from copy import deepcopy
import numpy as np

class MolecularCluster(AtomicCluster):
    def get_compare_exact(self, **kwargs):
        """this function quickly determines whether two clusters are identical
        given translational, rotational and permutational symmetries
        """
        permlist=self.get_permlist()
        return ExactMatchMolecularCluster(permlist=permlist, **kwargs)

    def get_mindist(self, **kwargs):
        """return a function which puts two structures in best alignment.
        
        take into account global rotational symmetry, global translational
        symmetry and permutational symmetry
        """
        permlist=self.get_permlist()
        return MinPermDistMolecularCluster(permlist=permlist, **kwargs)

class ExactMatchMolecularCluster(ExactMatchCluster):
    '''
    wrapper for ExactMatchCluster appropriate for molecules
    '''
    def __init__(self, permlist=None, can_invert=True, **kwargs):
        transform=TransformAtomicCluster(can_invert=can_invert)
        measure=MeasureMolecularCluster(permlist=permlist)

        ExactMatchCluster.__init__(self, transform=transform, measure=measure, **kwargs)

class MinPermDistMolecularCluster(MinPermDistCluster):
    '''
    wrapper for MinPermDistCluster appropriate for molecules
    '''
    def __init__(self, permlist=None, can_invert=True, **kwargs):
        transform=TransformAtomicCluster(can_invert=can_invert)
        measure=MeasureMolecularCluster(permlist=permlist)

        MinPermDistCluster.__init__(self, transform=transform, measure=measure, **kwargs)

class MeasureMolecularCluster(MeasureAtomicCluster):
    '''
    measure policies for molecular clusters; same as for atomic clusters but with
    different permutation rules
    '''
    def __init__(self, permlist=None):
        super(MeasureMolecularCluster, self).__init__(permlist)
    def find_permutation(self, X1, X2):
        '''
        overload atomic routine to properly permute molecules
        '''
        return find_best_permutation_molecular(X1, X2, self.permlist)

def find_best_permutation_molecular(X1, X2, permlist=None, user_algorithm=None,
                                        reshape=True, user_cost_matrix=_make_cost_matrix,
                                        **kwargs):
    '''
    modification of find_best_permutation that allows permutations of molecules,
    e.g. hydrogens are permuted with the oxygen to which they are bonded in a 
    water cluster
    '''
    if reshape:
        X1=X1.reshape([-1, 3])
        X2=X2.reshape([-1, 3])

    if permlist is None:
        permlist=[range(len(X1))]

    newperm=range(len(X1))
    disttot=0.

    for atomlists in permlist:
        if type(atomlists[0]) is not int:
            atomlist, associated=atomlists[0], atomlists[1:]
            sets=True
        else:
            atomlist=atomlists
            sets=False
        atomlist2=[newperm[a] for a in atomlist]
        if user_algorithm is not None:
            dist, perm=user_algorithm(X1[atomlist], X2[atomlist2], make_cost_matrix=user_cost_matrix, **kwargs)
        elif user_cost_matrix is not _make_cost_matrix:
            dist, perm=find_permutations_hungarian(X1[atomlist], X2[atomlist2], make_cost_matrix=user_cost_matrix, **kwargs)
        else:
            dist, perm=_find_permutations(X1[atomlist], X2[atomlist2], **kwargs)

        disttot += dist ** 2
        temp=deepcopy(newperm)
        for atom, i in zip(atomlist, xrange(len(atomlist))):
            temp[atom]=newperm[atomlist[perm[i]]]
        newperm=deepcopy(temp)
        if sets:
            for s in associated:
                temp=deepcopy(newperm)
                for atom, i in zip(s, xrange(len(s))):
                    temp[atom]=newperm[s[perm[i]]]
                newperm=deepcopy(temp)

    dist=np.sqrt(disttot)
    return dist, newperm

def permlist_water( nmol ):
    '''
    in this particular example, the coordinate array is expected to be ordered as follows:
        [ O1, O2, ..., ON, 1H1, 2H1, 1H2, 2H2, ..., 1HN, 2HN ] 
    '''
    permlist = [[range( nmol ), range( nmol, 3 * nmol, 2 ), range( nmol + 1, 3 * nmol, 2 )]]
    permlist += [list( x ) for x in zip( range( nmol, 3 * nmol, 2 ), range( nmol + 1, 3 * nmol, 2 ) )]
    return permlist