import numpy as np

from pele.systems import BaseSystem
from pele.systems import LJCluster
from pele.potentials._gupta import gupta

class GuptaCluster(LJCluster):
    def __init__(self, natoms, 
                 p=5.206, q=1.220, A=0.6124, xi=2.441):
        BaseSystem.__init__(self)

        self.potential_kwargs = dict(p=p,q=q,A=A,xi=xi)

        self.natoms = natoms

        self.params.database.accuracy = 1e-3
        self.params.basinhopping["temperature"] = 1.0
        self.params.gui.basinhopping_nsteps = 100

    def get_potential(self):
        return gupta(**self.potential_kwargs)


def rungui():
    from pele.gui import run_gui
    natoms = 13
    system = GuptaCluster(natoms)
    db = system.create_database()
    run_gui(system)

if __name__ == "__main__":
    rungui()
