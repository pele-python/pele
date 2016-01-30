from __future__ import division
import numpy as np
from pele.potentials import HS_WCA
from pele.potentials import InversePowerStillinger
from pele.optimize import LBFGS_CPP

class MinimizeUniformHardsSpheres(object):
    def __init__(self, nr_particles=42, hard_volume_fraction=0.5, epsilon=1, alpha=0.2, use_hswca=False):
        np.random.seed(42)
        self.nr_particles = nr_particles
        self.hard_volume_fraction = hard_volume_fraction
        self.epsilon = epsilon
        self.alpha = alpha
        self.use_hswca = use_hswca
        self.hard_radii = np.random.normal(loc=1, scale=0.1, size=self.nr_particles)
        self.box_length = np.power(np.sum(np.asarray([4 * np.pi * r**3 / 3 for r in self.hard_radii])) / self.hard_volume_fraction, 1/3)
        self.nr_dof = 3 * self.nr_particles
        self.x = np.random.uniform(-0.5 * self.box_length, 0.5 * self.box_length, self.nr_dof)
        self.box_vector = np.ones(3) * self.box_length
        self.rcut = 2 * (1 + alpha) * np.amax(self.hard_radii)
        if self.use_hswca:
            self.potential = HS_WCA(use_periodic=True, use_cell_lists=False, eps=self.epsilon, sca=self.alpha, radii=self.hard_radii, boxvec=self.box_vector, rcut=self.rcut)
        else:
            self.potential = InversePowerStillinger(8, boxvec=self.box_vector)
        self.optimizer = LBFGS_CPP(self.x, self.potential)
        print "energy before:", self.potential.getEnergy(self.x)
        self.optimizer.run()
        print "minimization converged", self.optimizer.get_result().success
        print "energy after:", self.potential.getEnergy(self.optimizer.get_result().coords)
        
if __name__ == "__main__":
    use_hswca = True
    for alpha in np.linspace(0, 0.5, 10):
        print "===, sca ==", alpha
        for phi in np.linspace(0.3, 0.6, 10):
            print "---"
            print "hard_volume_fraction:", phi
            print "soft_volume_fraction:", phi * (1 + alpha)**3
            print "( sca:", alpha, ")"
            MinimizeUniformHardsSpheres(hard_volume_fraction=phi, use_hswca=use_hswca)
