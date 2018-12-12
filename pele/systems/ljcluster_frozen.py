from __future__ import print_function
from __future__ import absolute_import
import numpy as np

from pele.systems import LJCluster
from pele.potentials import FrozenPotentialWrapper
from pele.mindist import optimize_permutations


__all__ = ["LJClusterFrozen"]


class LJClusterFrozen(LJCluster):
    def __init__(self, natoms, frozen_atoms, reference_coords):
        super(LJClusterFrozen, self).__init__(natoms)

        self.reference_coords = np.array(reference_coords)

        self.frozen_atoms = np.array(frozen_atoms)
        self.frozen_dof = np.array([list(range(3 * i, 3 * i + 3)) for i in self.frozen_atoms]).flatten()
        self.frozen_dof.sort()
        self.nfrozen = len(self.frozen_atoms)

        pot = self.get_potential()
        self.coords_converter = pot
        self.mobile_dof = self.coords_converter.get_mobile_dof()
        self.mobile_atoms = np.array([i for i in range(self.natoms) if i not in self.frozen_atoms], np.integer)
        self.nmobile = len(self.mobile_atoms)
        # print self.nmobile
        # print self.frozen_atoms
        # print self.mobile_dof


        if self.nfrozen <= 2:
            print("warning: Dealing properly with the rotational degrees of freedom and reflection symmetry in clusters with 1 or 2 frozen atoms is not implemented yet")
        if self.nfrozen == 3:
            print("warning: with three frozen atoms there is still reflection symmetry through the plane of the atoms.  This will not be accounted for.")


    def get_potential(self):
        pot = LJCluster.get_potential(self)
        frozen_pot = FrozenPotentialWrapper(pot, self.reference_coords, self.frozen_dof)
        return frozen_pot

    def get_permlist(self):
        """return the permutable mobile atoms"""
        # get permlist must be overloaded because the mindist functions will see the reduced set of coordinates
        return [list(range(self.nmobile))]

    def get_mindist(self):
        if self.get_nzero_modes() != 0:
            raise RuntimeError(
                "Dealing properly with the rotational degrees of freedom in clusters with 1 or 2 frozen atoms is not implemented yet")
        return lambda x1, x2: optimize_permutations(x1, x2, permlist=self.get_permlist())

    def get_orthogonalize_to_zero_eigenvectors(self):
        if self.get_nzero_modes() != 0:
            raise RuntimeError(
                "Dealing properly with the rotational degrees of freedom in clusters with 1 or 2 frozen atoms is not implemented yet")
        return None

    # return lambda grad, coords: grad
    # return SubtractZeroEigs(self.data.frozenlist)

    def get_compare_exact(self):
        mindist = self.get_mindist()
        accuracy = 0.01
        return lambda x1, x2: mindist(x1, x2)[0] < accuracy

    def get_pgorder(self, coords):
        if self.get_nzero_modes() != 0:
            raise RuntimeError(
                "Dealing properly with the rotational degrees of freedom in clusters with 1 or 2 frozen atoms is not implemented yet")
        return 1

    def get_random_configuration(self):
        """return a random configuration in the reduced coordinates"""
        boxl = 0.7 * float(self.nmobile) ** (1. / 3)
        coords = np.random.uniform(-1, 1, [3 * self.nmobile]) * boxl
        return coords

    def get_nzero_modes(self):
        if self.nfrozen == 0:
            # translational + rotational degrees of freedom
            return 6
        if self.nfrozen == 1:
            # only rotational
            return 3
        if self.nfrozen == 2:
            # still have 1 free rotational degree of freedom 
            return 1
        return 0


    def draw(self, coordslinear, index):
        """draw the frozen atoms differently from the mobile atoms"""
        from ._opengl_tools import draw_atomic_binary

        full_coords = self.coords_converter.get_full_coords(coordslinear)
        draw_atomic_binary(full_coords, index, self.mobile_atoms,
                           self.frozen_atoms, subtract_com=False, rA=0.5,
                           rB=0.5)


#
# testing only below here
#

def test():  # pragma: no cover
    from pele.gui import run_gui

    natoms = 13

    fsys = LJCluster(natoms)
    reference_coords = fsys.get_random_configuration()
    frozen_atoms = [0, 2, 3, 4]
    system = LJClusterFrozen(natoms, frozen_atoms, reference_coords)

    run_gui(system)


if __name__ == "__main__":
    test()

