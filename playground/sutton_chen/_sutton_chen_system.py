import numpy as np

from pele.systems import BaseSystem
from pele.systems.morse_bulk import MorseBulk
from pele.potentials._sutton_chen import SuttonChen

class SuttonChenBulk(MorseBulk):
    def __init__(self, natoms, eps=1., sig=1., c=1., boxvec=[10., 10., 10.], rcut=2.5,
                  n=10, m=8):
        BaseSystem.__init__(self)
        
        self.potential_kwargs = dict(eps=eps, sig=sig, c=c, boxvec=boxvec, rcut=rcut, n=n, m=m)
        
        self.natoms = natoms
        self.boxvec = np.array(boxvec, dtype=float)
        self.periodic = True
        
        self.r0 = sig # the equilibrium separation of the atoms.

        self.params.database.accuracy = 1e-3
        self.params.basinhopping["temperature"] = 1.0
        self.params.gui.basinhopping_nsteps = 100
    
    def get_potential(self):
        return SuttonChen(**self.potential_kwargs)


def rungui():
    from pele.gui import run_gui
    natoms = 17
    boxl = 10.
    boxvec = np.ones(3) * boxl
    system = SuttonChenBulk(natoms, rcut=30000.2, boxvec=boxvec, c=144.41, n=12, m=6)
#    system = MorseBulk(natoms, boxvec, rho=3., r0=1., A=1.)
    db = system.create_database()
    run_gui(system)

if __name__ == "__main__":
    rungui()
