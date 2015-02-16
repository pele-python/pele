import numpy as np

from pele.systems import BaseSystem
from pele.mindist.periodic_exact_match import ExactMatchPeriodic, MeasurePeriodic
from pele.mindist.periodic_mindist import MinDistBulk
from pele.systems.morse_bulk import put_in_box
from pele.potentials import InversePower


class SoftSpereSystem(BaseSystem):
    """Binary Lennard Jones potential with periodic boundary conditions"""

    def __init__(self, radii, boxvec, power=2):
        super(SoftSpereSystem, self).__init__(self)

        self.radii = np.array(radii)
        self.natoms = self.radii.size 
        self.boxvec = np.array(boxvec)
        self.periodic = True
        self.potential_kwargs = dict()
        self.power = power
        self.eps = 1.
        
        self.set_params(self.params)

    def set_params(self, params):
        params.database.accuracy = 1e-4
        params.structural_quench_params.tol = 1e-6
        
        params.double_ended_connect.local_connect_params.tsSearchParams.iprint = 1
        params.double_ended_connect.local_connect_params.tsSearchParams.hessian_diagonalization = True

    def get_potential(self):
        return InversePower(self.power, self.eps, self.radii, boxvec=self.boxvec, **self.potential_kwargs)

    def get_random_configuration(self):
        x = np.zeros([self.natoms, self.boxvec.size])
        for i in range(self.boxvec.size):
            x[:, i] = np.random.uniform(-self.boxvec[i] / 2., self.boxvec[i] / 2., self.natoms)
        return x.ravel()

    def draw(self, coordslinear, index):
        from pele.systems._opengl_tools import draw_atomic_binary_polydisperse, draw_box
        put_in_box(coordslinear, self.boxvec)
        draw_atomic_binary_polydisperse(coordslinear, index, bdim=self.boxvec.size, radii=self.radii)
        draw_box(self.boxvec)

    def get_permlist(self):
        return []

    def get_mindist(self):
        measure = MeasurePeriodic(self.boxvec, permlist=self.get_permlist())
        return MinDistBulk(self.boxvec, measure)
    
    def get_compare_exact(self):
        accuracy = self.params.database.accuracy
        measure = MeasurePeriodic(self.boxvec, self.get_permlist())
        compare = ExactMatchPeriodic(measure, accuracy=accuracy)
        return compare

    def get_orthogonalize_to_zero_eigenvectors(self):
        # TODO: there are some zero eigenvectors which can be removed 
        return None
    
    def get_system_properties(self):
        return dict(natoms=int(self.natoms),
                    radii=self.radii,
                    power=self.power,
                    eps=self.eps,
                    boxvec=self.boxvec,
                    potential="Soft Sphere Bulk",
                    potential_kwargs=self.potential_kwargs,
        )
    
    def get_metric_tensor(self, coords):
        return None



def rungui():  # pragma: no cover
    from pele.gui import run_gui

    natoms = 17
    boxl = 2.
    boxvec = np.ones(3) * boxl
    # system = MorseCluster(natoms, rho=1.6047, r0=2.8970, A=0.7102, rcut=9.5)
    radii = np.ones(natoms) * .8
    system = SoftSpereSystem(radii, boxvec)
    db = system.create_database()
    run_gui(system, db)


if __name__ == "__main__":
    rungui()
