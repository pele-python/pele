import numpy as np
from copy import copy
from pygmin.potentials import BasePotential

import networkx as nx

__all__ = ["XYModel"]

class RectangularLattice(object):
    """
    for getting connectivity of a rectangular lattice
    """
    def __init__(self, Lx, Ly):
        self.Lx = Lx
        self.Ly = Ly
        self.nspins = self.Lx * self.Ly
                            
    def i2xy(self, i):
        #xy = np.zeros(2)
        i -= self.nspins * int(np.floor( float(i) / self.nspins ))
        
        x = i % self.Lx
        y = int(np.floor(i / self.Lx))
        return x, y
    
    def xy2i(self, xy):
        #xy = np.zeros(2)
        x = xy[0]
        y = xy[1]
        x -= self.Lx * int(np.floor( float(x) / self.Lx ))
        y -= self.Ly * int(np.floor( float(y) / self.Ly ))
        
        i = x + y * self.Lx
        return i
        
    
    

class XYModel(BasePotential):
    """
    XY model of 2d spins on a lattice
    """
    def __init__(self, dim = [4, 4], phi=np.pi):
        self.dim = copy(dim)
        self.nspins = np.prod(dim)
        
        self.G = nx.grid_graph(dim, periodic=True)
        
        self.phases = dict()
        binary_disorder = True
        if binary_disorder:
            for edge in self.G.edges():
                self.phases[edge] = phi * np.random.random_integers(0,1)
        else:
            for edge in self.G.edges():
                self.phases[edge] = np.random.uniform(-phi, phi)
        nx.set_edge_attributes(self.G, "phase", self.phases)

        self.indices = dict()
        i = 0
        for node in self.G.nodes():
            self.indices[node] = i
            i += 1 
        
        self.num_edges = self.G.number_of_edges()

        
        
    def getEnergy(self, angles):
        #do internal energies first
        E = 0.
        for edge in self.G.edges():
            phase = self.phases[edge]
            u = self.indices[edge[0]]
            v = self.indices[edge[1]]
            E += np.cos( -angles[u] + angles[v] + phase )
        #E = self.num_edges - E
        E = - E
        return E
        
    def getEnergyGradient(self, angles):
        #do internal energies first
        E = 0.
        grad = np.zeros(self.nspins)
        for edge in self.G.edges():
            phase = self.phases[edge]
            u = self.indices[edge[0]]
            v = self.indices[edge[1]]
            E += np.cos( -angles[u] + angles[v] + phase )
            
            g = -np.sin( -angles[u] + angles[v] + phase )
            grad[u] += g
            grad[v] += -g
        #E = self.num_edges - E
        E =  - E
        return E, grad



def test_basin_hopping(pot, angles):
    from pygmin.basinhopping import BasinHopping
    from pygmin.takestep.displace import RandomDisplacement
    from pygmin.takestep.adaptive import AdaptiveStepsize
    
    takestep = RandomDisplacement(stepsize = np.pi/4)
    takestepa = AdaptiveStepsize(takestep, frequency = 20)
    
    bh = BasinHopping( angles, pot, takestepa, temperature = 1.01)
    bh.run(20)

def test():
    pi = np.pi
    L = 3
    nspins = L**2
    
    #phases = np.zeros(nspins)
    pot = XYModel( dim = [L,L], phi = np.pi) #, phases=phases)
    
    
    angles = np.random.uniform(-pi, pi, nspins)
    print angles

    e = pot.getEnergy(angles)
    print "energy ", e

    print "numerical gradient"
    ret = pot.getEnergyGradientNumerical(angles)
    print ret[1]
    print "analytical gradient"
    ret2 = pot.getEnergyGradient(angles)
    print ret2[1]
    print ret[0]
    print ret2[0]
    

    
    #try a quench
    from pygmin.optimize import mylbfgs
    ret = mylbfgs(angles, pot)
    
    print "quenched e = ", ret.energy
    print ret.coords
    
    test_basin_hopping(pot, angles)

if __name__ == "__main__":
    test()
