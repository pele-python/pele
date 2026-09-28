from itertools import groupby

import numpy as np

from pele.optimize import modifiedfire_cpp
from pele.potentials import LAMMPSPotential
from pele.systems import AtomicCluster, dict_copy_update
from pele.utils.elements import elements


class AtomicClusterLAMMPS(AtomicCluster):
    def __init__(self, lmp):
        super().__init__()
        self.lmp = lmp
        self.natoms = lmp.get_natoms()
        self.params.database.accuracy = 1e-5
        self.params.basinhopping["temperature"] = 1.0

    def get_permlist(self):
        atom_types = self.lmp.extract_atom("type")
        atom_numbers = sorted(range(self.natoms), key=atom_types.__getitem__)
        return [list(group) for key, group in groupby(atom_numbers, key=atom_types.__getitem__)]

    def get_potential(self):
        return LAMMPSPotential(self.lmp)

    def get_system_properties(self):
        return dict(natoms=int(self.natoms),
                    potential="LAMMPS cluster",
                    )

    def get_minimizer(self, **kwargs):
        pot = self.get_potential()
        kwargs = dict_copy_update(self.params["structural_quench_params"], kwargs)
        return lambda coords: modifiedfire_cpp(coords, pot, **kwargs)

    def draw(self, coordslinear, index, subtract_com=True):
        from pele.systems._opengl_tools import draw_sphere, change_color
        coords = coordslinear.reshape(-1, 3)
        if subtract_com:
            com = np.mean(coords, axis=0)
            coords = coords - com[np.newaxis, :]
        for x, t in zip(coords, self.lmp.numpy.extract_atom("type")):
            color = elements[t + 3]["color"]
            change_color(color)
            draw_sphere(x, radius=0.5)
