from __future__ import absolute_import
import numpy as np
from copy import copy
from pele.potentials import BasePotential

import networkx as nx

__all__ = ["XYModel"]


def angle_to_2dvector(theta):
    return np.cos(theta), np.sin(theta)


class XYModel(BasePotential):
    """
    XY model of 2d spins on a lattice
    """

    def __init__(self, dim=None, phi=np.pi, periodic=True, phases=None):
        if not dim: dim = [4, 4]
        dim = copy(dim)
        self.dim = copy(dim)
        self.nspins = np.prod(dim)

        self.G = nx.grid_graph(dim, periodic)

        if phases is not None:
            self.phases = phases
        else:
            self.phases = dict()
            binary_disorder = True
            if binary_disorder:
                for edge in self.G.edges():
                    self.phases[edge] = phi * np.random.random_integers(0, 1)
            else:
                for edge in self.G.edges():
                    self.phases[edge] = np.random.uniform(-phi, phi)
        nx.set_edge_attributes(self.G, self.phases, "phase")

        self.indices = dict()
        self.index2node = dict()
        nodes = sorted(self.G.nodes())
        for i, node in enumerate(nodes):
            self.indices[node] = i
            self.index2node[i] = node

        self.num_edges = self.G.number_of_edges()

        self.set_up_neighborlists()

    def get_phases(self):
        return self.phases.copy()

    def set_up_neighborlists(self):
        neighbors = []
        self.phase_matrix = np.zeros([self.nspins, self.nspins])
        for edge in self.G.edges():
            u = self.indices[edge[0]]
            v = self.indices[edge[1]]
            neighbors.append([u, v])
            self.phase_matrix[u, v] = self.phases[edge]
            self.phase_matrix[v, u] = self.phases[edge]

        self.neighbors = np.array(neighbors).reshape([-1, 2])

    def get_spin_energies(self, angles):
        """return the local energy of each spin"""
        energies = np.zeros(angles.size)
        for edge in self.G.edges():
            phase = self.phases[edge]
            u = self.indices[edge[0]]
            v = self.indices[edge[1]]
            E = -np.cos(-angles[u] + angles[v] + phase)
            energies[u] += E
            energies[v] += E
        return energies

    def getEnergy(self, angles):
        e, g = self.getEnergyGradient(angles)
        return e

    def getEnergyGradient(self, angles):
        from . import _cython_tools

        return _cython_tools.xymodel_energy_gradient(angles, self.phase_matrix, self.neighbors)

    # def getEnergyGradient(self, angles):

# # do internal energies first
# E = 0.
#        grad = np.zeros(self.nspins)
#        for edge in self.G.edges():
#            phase = self.phases[edge]
#            u = self.indices[edge[0]]
#            v = self.indices[edge[1]]
#            E += np.cos( -angles[u] + angles[v] + phase )
#            
#            g = -np.sin( -angles[u] + angles[v] + phase )
#            grad[u] += g
#            grad[v] += -g
#        E =  - E
#        return E, grad



#def test_basin_hopping(pot, angles):
#    from pele.basinhopping import BasinHopping
#    from pele.takestep.displace import RandomDisplacement
#    from pele.takestep.adaptive import AdaptiveStepsize
#    
#    takestep = RandomDisplacement(stepsize = np.pi/4)
#    takestepa = AdaptiveStepsize(takestep, frequency = 20)
#    
#    bh = BasinHopping( angles, pot, takestepa, temperature = 1.01)
#    bh.run(20)
#
#def test():
#    pi = np.pi
#    L = 3
#    nspins = L**2
#    
#    #phases = np.zeros(nspins)
#    pot = XYModel( dim = [L,L], phi = np.pi) #, phases=phases)
#    
#    
#    angles = np.random.uniform(-pi, pi, nspins)
#    print angles
#
#    e = pot.getEnergy(angles)
#    print "energy ", e
#
#    print "numerical gradient"
#    ret = pot.getEnergyGradientNumerical(angles)
#    print ret[1]
#    print "analytical gradient"
#    ret2 = pot.getEnergyGradient(angles)
#    print ret2[1]
#    print ret[0]
#    print ret2[0]
#    
#
#    
#    #try a quench
#    from pele.optimize import mylbfgs
#    ret = mylbfgs(angles, pot)
#    
#    print "quenched e = ", ret.energy
#    print ret.coords
#    
#    test_basin_hopping(pot, angles)
#
#if __name__ == "__main__":
#    test()

