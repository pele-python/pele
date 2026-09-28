import numpy as np
from lammps import lammps, LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR

from pele.potentials import BasePotential


class LAMMPSPotential(BasePotential):
    def __init__(self, lmp: lammps) -> None:
        self.lmp = lmp
        lmp.cmd.neigh_modify(every=1, delay=0, check=True)
        self.coords = lmp.numpy.extract_atom('x')
        self.forces = lmp.numpy.extract_atom('f')

    def getEnergy(self, coords: np.ndarray) -> float:
        self.coords[:] = coords.reshape(-1, 3)
        self.lmp.cmd.run(0)
        return self.lmp.extract_compute("thermo_pe", LMP_STYLE_GLOBAL, LMP_TYPE_SCALAR)

    def getEnergyGradient(self, coords: np.ndarray) -> (float, np.ndarray):
        energy = self.getEnergy(coords)
        return energy, -self.forces.flatten()
