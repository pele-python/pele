from __future__ import print_function
from __future__ import absolute_import
import numpy as np

from pele.systems import BLJBulk
from pele.potentials import FrozenPotentialWrapper
from pele.mindist import optimize_permutations
from pele.systems.morse_bulk import put_in_box


class BLJBulkFrozen(BLJBulk):
    """Binary Lennard Jones in a periodic box with frozen atoms"""
    def __init__(self, natoms, boxvec, reference_coords, frozen_atoms, 
                  ntypeA="default", **potential_kwargs):
        super(BLJBulkFrozen, self).__init__(natoms, boxvec, ntypeA=ntypeA, **potential_kwargs)
        ntypeA = self.ntypeA
        
        self.reference_coords = np.array(reference_coords)

        self.frozen_atoms = np.array(frozen_atoms, dtype=int)
        self.frozen_dof = np.array([list(range(3 * i, 3 * i + 3)) for i in self.frozen_atoms]).flatten()
        self.frozen_dof.sort()
        self.nfrozen = len(self.frozen_atoms)

        pot = self.get_potential()
        self.coords_converter = pot
        self.mobile_dof = self.coords_converter.get_mobile_dof()
        self.mobile_atoms = np.array([i for i in range(self.natoms) if i not in self.frozen_atoms], np.integer)
        self.nmobile = len(self.mobile_atoms)
        
        self.mobile_Aatoms = [i for i in self.mobile_atoms if i < ntypeA]
        self.mobile_Batoms = [i for i in self.mobile_atoms if i >= ntypeA]
        self.frozen_Aatoms = [i for i in self.frozen_atoms if i < ntypeA]
        self.frozen_Batoms = [i for i in self.frozen_atoms if i >= ntypeA]

    def get_potential(self):
        try:
            return self.potential
        except AttributeError:
            pass
        blj = super(BLJBulkFrozen, self).get_potential()
        self.potential = FrozenPotentialWrapper(blj, self.reference_coords, self.frozen_dof)
        return self.potential

    def get_permlist(self):
        """return the permutable mobile atoms"""
        # get permlist must be overloaded because the mindist functions will see the reduced set of coordinates
        return [list(range(len(self.mobile_Aatoms))), list(range(len(self.mobile_Batoms)))]

    def get_mindist(self):
        return lambda x1, x2: optimize_permutations(x1, x2, permlist=self.get_permlist())

    def get_orthogonalize_to_zero_eigenvectors(self):
        return None

    def get_compare_exact(self):
        mindist = self.get_mindist()
        accuracy = 0.01
        return lambda x1, x2: mindist(x1, x2)[0] < accuracy

    def get_pgorder(self, coords):
        return 1

    def get_random_configuration(self):
        """return a random configuration in the reduced coordinates"""
        boxl = 0.7 * float(self.nmobile) ** (1. / 3)
        coords = np.random.uniform(-1, 1, [3 * self.nmobile]) * boxl
        return coords

    def get_nzero_modes(self):
        return 0

    def get_system_properties(self):
        return dict(natoms=int(self.natoms),
                    ntypeA=int(self.ntypeA),
                    boxvec=self.boxvec,
                    potential="BLJ Bulk Frozen",
                    potential_kwargs=self.potential_kwargs,
                    reference_coords=self.reference_coords,
                    frozen_atoms=self.frozen_atoms,
        )

    def draw(self, coordslinear, index):
        """draw the frozen atoms differently from the mobile atoms"""
        from ._opengl_tools import draw_atoms, draw_box

        full_coords = self.coords_converter.get_full_coords(coordslinear)
        put_in_box(full_coords, self.boxvec)
        
        rA = 0.5
        rB = 0.44
        
        # draw A atoms
        if index == 1:
            cA = [0.65, 0.0, 0.0, 1.]
            cB = [0.25, 0.00, 0., 1.]
            cF = [0.24, 0.25, 0.25, 1.]
        else:
            cA = cB = cF = [0.00, 0.65, 0., 1.]
        draw_atoms(full_coords, self.mobile_Aatoms, radius=rA, color=cA)
        draw_atoms(full_coords, self.frozen_Aatoms, radius=rA, color=cF)
        draw_atoms(full_coords, self.mobile_Batoms, radius=rB, color=cB)
        draw_atoms(full_coords, self.frozen_Batoms, radius=rB, color=cF)
        
        draw_box(self.boxvec)


#
# testing only below here
#

def rungui():  # pragma: no cover
    from pele.gui import run_gui

    natoms = 40
    boxl = 3.
    boxvec = np.ones(3) * boxl
    fsys = BLJBulk(natoms, boxvec)
    # system = MorseCluster(natoms, rho=1.6047, r0=2.8970, A=0.7102, rcut=9.5)
    reference_coords = fsys.get_random_configuration()
    frozen_atoms = [0, 2, 3, 4, natoms-1]

    system = BLJBulkFrozen(natoms, boxvec, reference_coords, frozen_atoms)
    print(system.ntypeA)
    db = system.create_database()
    run_gui(system, db)


if __name__ == "__main__":
    rungui()

