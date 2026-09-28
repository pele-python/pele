import numpy as np
from lammps import lammps, LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR

from pele.potentials import BasePotential


def prepare_lammps(lmp: lammps) -> None:
    """set up a LAMMPS instance for use as a pele potential

    pele passes coordinates and gradients in LAMMPS's local atom order, so
    that order must not change: turn off LAMMPS's spatial sorting of atoms,
    and require a single MPI rank (with more, each rank holds only its own
    atoms). Rebuild neighbor lists whenever needed, since a minimizer can move
    atoms arbitrarily far between calls.
    """
    if lmp.extract_setting("world_size") > 1:
        raise ValueError("LAMMPSPotential needs LAMMPS running on a single MPI rank")
    lmp.command("atom_modify sort 0 0.0")
    lmp.command("neigh_modify every 1 delay 0 check yes")


class LAMMPSPotential(BasePotential):
    def __init__(self, lmp: lammps) -> None:
        prepare_lammps(lmp)
        self.lmp = lmp
        self.coords = lmp.numpy.extract_atom('x')
        self.forces = lmp.numpy.extract_atom('f')

    def getEnergy(self, coords: np.ndarray) -> float:
        self.coords[:] = coords.reshape(-1, 3)
        self.lmp.cmd.run(0)
        return self.lmp.extract_compute("thermo_pe", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR)

    def getEnergyGradient(self, coords: np.ndarray) -> (float, np.ndarray):
        energy = self.getEnergy(coords)
        return energy, -self.forces.flatten()
