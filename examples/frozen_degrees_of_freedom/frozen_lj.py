"""
this example shows how to freeze degrees of freedom using the Lennard Jones potential as
an example
"""
import numpy as np
from pele.potentials import LJ, FrozenPotentialWrapper
from pele.optimize import mylbfgs


def main():
    natoms = 4
    pot = LJ()

    reference_coords = np.random.uniform(-1, 1, [3 * natoms])
    print reference_coords

    # freeze the first two atoms (6 degrees of freedom)
    frozen_dof = range(6)

    fpot = FrozenPotentialWrapper(pot, reference_coords, frozen_dof)

    reduced_coords = fpot.get_reduced_coords(reference_coords)

    print "the energy in the full representation:"
    print pot.getEnergy(reference_coords)
    print "is the same as the energy in the reduced representation:"
    print fpot.getEnergy(reduced_coords)

    ret = mylbfgs(reduced_coords, fpot)
    print "after a minimization the energy is ", ret.energy, "and the rms gradient is", ret.rms
    print "the coordinates of the frozen degrees of freedom are unchanged"
    print "starting coords:", reference_coords
    print "minimized coords:", fpot.get_full_coords(ret.coords)


if __name__ == "__main__":
    main()