import numpy as np

from pele.systems import BaseSystem
from pele.systems import LJCluster
from pele.potentials._fin_sin import FinSin

class FinSinCluster(LJCluster):
    def __init__(self, natoms, A=2.010637):
        BaseSystem.__init__(self)
        
        self.potential_kwargs = dict(A=A)
        
        self.natoms = natoms
        
        self.params.database.accuracy = 1e-3
        self.params.basinhopping["temperature"] = 1.0
        self.params.gui.basinhopping_nsteps = 100
    
    def get_potential(self):
        return FinSin(**self.potential_kwargs)


def rungui():
    from pele.gui import run_gui
    natoms = 17
    system = FinSinCluster(natoms)
    db = system.create_database()
    run_gui(system)

if __name__ == "__main__":
    rungui()
