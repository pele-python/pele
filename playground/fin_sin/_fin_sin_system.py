import numpy as np

from pele.systems import BaseSystem
from pele.systems import LJCluster
from pele.potentials._fin_sin import FinSin

class FinSinCluster(LJCluster):
    def __init__(self, natoms, d=4.114825, A=1.887117, beta=0.0, c=3.25,
                 c0=43.4475218, c1=-31.9332978, c2=6.0804249):
        BaseSystem.__init__(self)
        
        self.potential_kwargs = dict(d=d, A=A, beta=beta, c=c, c0=c0, c1=c1, c2=c2)
        
        self.natoms = natoms
        
        self.params.database.accuracy = 1e-3
        self.params.basinhopping["temperature"] = 1.0
        self.params.gui.basinhopping_nsteps = 100
    
    def get_potential(self):
        return FinSin(**self.potential_kwargs)

    def get_system_properties(self):
        return dict(natoms=int(self.natoms),
                    potential="Finnis-Sinclair cluster",
                    potential_kwargs=self.potential_kwargs
                    )

    def draw(self, coordslinear, index): # pragma: no cover
        """
        tell the gui how to represent your system using openGL objects
        
        Parameters
        ----------
        coords : array
        index : int
            we can have more than one molecule on the screen at one time.  index tells
            which one to draw.  They are viewed at the same time, so they should be
            visually distinct, e.g. different colors.  accepted values are 1 or 2        
        """
        from pele.systems._opengl_tools import draw_atomic_single_atomtype
#        radius = self.potential_kwargs["d"] / 4
        radius = 1.
        draw_atomic_single_atomtype(coordslinear, index, subtract_com=True,
                                    radius=radius)


def rungui():
    from pele.gui import run_gui
    natoms = 13
    system = FinSinCluster(natoms)
    db = system.create_database()
    run_gui(system)

if __name__ == "__main__":
    rungui()
