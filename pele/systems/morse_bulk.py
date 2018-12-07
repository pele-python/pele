from __future__ import absolute_import
import numpy as np

from .morse_cluster import MorseCluster
from pele.potentials import Morse
from pele.mindist.periodic_exact_match import ExactMatchPeriodic, MeasurePeriodic
from pele.mindist import optimize_permutations
from pele.transition_states import InterpolateLinearMeasure

def put_in_box(x, boxvec):
    x = x.reshape(-1, boxvec.size)
    x -= boxvec * np.round(x / boxvec)


class MorseBulk(MorseCluster):
    """morse potential with periodic boundary conditions"""

    def __init__(self, natoms, boxvec, rho=2., r0=1., A=1., rcut=None):
        super(MorseBulk, self).__init__(natoms, rho=rho, r0=r0, A=A, rcut=rcut)

        self.boxvec = boxvec
        self.periodic = True
        
        self.params.double_ended_connect.local_connect_params.NEBparams.interpolator = InterpolateLinearMeasure(MeasurePeriodic(self.boxvec))


    def get_random_configuration(self):
        x = np.zeros([self.natoms, 3])
        for i in range(3):
            x[:, i] = np.random.uniform(-self.boxvec[i] / 2., self.boxvec[i] / 2., self.natoms)
        return x.flatten()

    def get_potential(self):
        return Morse(rho=self.rho, r0=self.r0, A=self.A, boxvec=self.boxvec, rcut=self.rcut)

    def draw(self, coordslinear, index):
        from pele.systems._opengl_tools import draw_box
        put_in_box(coordslinear, self.boxvec)
        draw_box(self.boxvec)
        MorseCluster.draw(self, coordslinear, index, subtract_com=False)

    def get_mindist(self):
        return lambda x1, x2: optimize_permutations(x1, x2, permlist=self.get_permlist(),
                                                    box_lengths=self.boxvec)

    def get_compare_exact(self):
        accuracy = self.params.database.accuracy
        measure = MeasurePeriodic(self.boxvec, self.get_permlist())
        compare = ExactMatchPeriodic(measure, accuracy=accuracy)
        return compare

    def get_orthogonalize_to_zero_eigenvectors(self):
        # TODO: there are some zero eigenvectors which can be removed 
        return None


def rungui():  # pragma: no cover
    from pele.gui import run_gui

    natoms = 17
    boxl = 2.
    boxvec = np.ones(3) * boxl
    # system = MorseCluster(natoms, rho=1.6047, r0=2.8970, A=0.7102, rcut=9.5)
    system = MorseBulk(natoms, boxvec, rho=3., r0=1., A=1.)
    db = system.create_database()
    run_gui(system, db)


if __name__ == "__main__":
    rungui()

