"""
this example shows how to freeze degrees of freedom using the Lennard Jones potential as
an example
"""
import numpy as np
from pele.potentials import LJ
from pele.utils.frozen_atoms import FrozenPotWrapper
from pele.optimize import mylbfgs
from pele.systems import LJCluster
from pele.mindist import optimize_permutations

class LJClusterFrozen2D(LJCluster):
    def __init__(self, natoms):
        super(LJClusterFrozen2D, self).__init__(natoms)
        
        self.frozen_dof = np.array(range(2,3*natoms,3))
        assert len(self.frozen_dof) == self.natoms
        print "frozen dof", self.frozen_dof
        
        x = self.get_random_configuration()
        self.reference_coords = np.zeros(3*self.natoms)
        self.reference_coords[0::3] = x[0::2]
        self.reference_coords[1::3] = x[1::2]
        print "reference coords", self.reference_coords
        
        pot = self.get_potential()
        self.coords_converter = pot.coords_converter
        self.mobile_dof = self.coords_converter.get_mobile_dof()
        print "mobile dof", self.mobile_dof
#        print self.nmobile
#        print self.frozen_atoms
#        print self.mobile_dof


        if len(self.frozen_dof) <= 1:
            print "warning: Dealing properly with the rotational and translational degrees of freedom and reflection symmetry in clusters with 1 or 2 frozen atoms is not implemented yet"
            
    def get_potential(self):
        pot = LJCluster.get_potential(self)
        frozen_pot = FrozenPotWrapper(pot, self.reference_coords, self.frozen_dof)
        return frozen_pot
    
    def get_mindist(self):
        def mindist2d(x1, x2):
            mdist3d = super(LJClusterFrozen2D, self).get_mindist()
            x1full = self.coords_converter.get_full_coords(x1)
            x2full = self.coords_converter.get_full_coords(x2)
            dist, x1new_full, x2new_full = mdist3d(x1full, x2full)
            x1new = self.coords_converter.get_reduced_coords(x1new_full)
            x2new = self.coords_converter.get_reduced_coords(x2new_full)
            return dist, x1new, x2new
            
        return mindist2d

    def get_orthogonalize_to_zero_eigenvectors(self):
#         if self.get_nzero_modes() != 0:
#             raise RuntimeError("Dealing properly with the rotational degrees of freedom in clusters with 1 or 2 frozen atoms is not implemented yet") 
        return None
#        return lambda grad, coords: grad
#        return SubtractZeroEigs(self.data.frozenlist)

    def get_compare_exact(self):
        mindist = self.get_mindist()
        accuracy = 0.01
        return lambda x1, x2: mindist(x1, x2)[0] < accuracy
 
    def get_pgorder(self, coords):
        print "warning point group order is not correct for the 2d system"
        return 1

    def get_random_configuration(self):
        """return a random configuration in the reduced coordinates"""
        boxl = 0.7 * float(self.natoms)**(1./2)
        coords = np.random.uniform(-1, 1, [2*self.natoms]) * boxl
        return coords

    def get_nzero_modes(self):
        # 2 translational + 1 rotational
        return 3



    def draw(self, coordslinear, index):
        """draw the frozen atoms differently from the mobile atoms"""
        from pele.systems._opengl_tools import draw_atomic_single_atomtype
        full_coords = self.coords_converter.get_full_coords(coordslinear)
        draw_atomic_single_atomtype(full_coords, index, subtract_com=False)


def main():
    natoms = 8
    system = LJClusterFrozen2D(natoms)
    
    db = system.create_database()
    
    from pele.gui import run_gui
    run_gui(system, db=db)
    
    

if __name__ == "__main__":
    main()