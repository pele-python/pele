from __future__ import print_function
import numpy as np
from copy import copy
from numpy import sin, cos
import networkx as nx

from pele.potentials import BasePotential
import pele.utils.rotations as rotations

__all__ = ["HeisenbergModel"]


def make3dVector(u):
    """
    make a 3d unit vector from (theta, phi)
    """
    sinphi = sin(u[1])

    vec = np.zeros(3)
    vec[0] = sinphi * cos(u[0])
    vec[1] = sinphi * sin(u[0])
    vec[2] = cos(u[1])
    if np.abs(np.linalg.norm(vec) - 1) > 1e-5:
        print("make3dVector: vector not normalized", u, vec, np.linalg.norm(vec))
    return vec


def make2dVector(u):
    """
    make (theta, phi) from a 3d vector
    """
    vec = np.zeros(2)
    vec[1] = np.arccos(u[2])
    vec[0] = np.arctan2(u[1], u[0])
    return vec


def coords2ToCoords3(coords2):
    if len(np.shape(coords2)) == 1:
        nvec = len(coords2) // 2
        coords2 = np.reshape(coords2, [nvec, 2])
    else:
        nvec = len(coords2[:, 0])
    coords3 = np.zeros([nvec, 3])
    for i in range(nvec):
        coords3[i, :] = make3dVector(coords2[i, :])
    return coords3


def coords3ToCoords2(coords3):
    if len(np.shape(coords3)) == 1:
        nvec = len(coords3) / 3
        coords3 = np.reshape(coords3, [nvec, 3])
    else:
        nvec = len(coords3[:, 0])
    coords2 = np.zeros([nvec, 2])
    for i in range(nvec):
        coords2[i, :] = make2dVector(coords3[i, :])
    return coords2


def makeGrad2(vec2, grad3):
    grad2 = np.zeros(2)
    c0 = cos(vec2[0])
    c1 = cos(vec2[1])
    s0 = sin(vec2[0])
    s1 = sin(vec2[1])
    grad2[0] = -s0 * grad3[0] + c0 * grad3[1]
    grad2[1] = c0 * c1 * grad3[0] + s0 * c1 * grad3[1] - s1 * grad3[2]
    grad2[0] *= s1  # I need this to agree with the numerical gradient, but I think it shouldn't be there
    return grad2


def grad3ToGrad2(coords2, grad3):
    if len(np.shape(grad3)) == 1:
        nvec = len(grad3) / 3
        grad3 = np.reshape(grad3, [nvec, 3])
    else:
        nvec = len(grad3[:, 0])
    if len(np.shape(coords2)) == 1:
        coords2 = np.reshape(coords2, [nvec, 2])
    grad2 = np.zeros([nvec, 2])
    for i in range(nvec):
        grad2[i, :] = makeGrad2(coords2[i, :], grad3[i, :])
    return grad2


class HeisenbergModel(BasePotential):
    """
    The classical Heisenberg Model of 3d spins on a lattice.
    
    Parameters
    ----------
    dim: list
        an array giving the dimensions of the lattice
    field_disorder: float
        the magnitude of the randomness in the fields
    fields: array
        use this to explicitly set the values of the field disorder
    
    Notes
    -----
    The Hamiltonian is::
    
        H = - sum_ij J dot( s_i, s_j )
    
    where s_i are normalized 3d vectors.
    
    This can be generalized for disordered systems to::
    
        H = - sum_ij J_ij dot( s_i, s_j )  - sum_i dot( h_i, s_i )
    
    where h_i are quenched random variables.  (h_i is a vector)
    """

    def __init__(self, dim=None, field_disorder=1., fields=None):
        if dim is None: dim = [4, 4]
        self.dim = copy(dim)
        self.nspins = np.prod(dim)

        self.G = nx.grid_graph(dim, periodic=True)

        self.fields = np.zeros([self.nspins, 3])

        self.indices = dict()
        nodes = sorted(self.G.nodes())
        for i, node in enumerate(nodes):
            self.indices[node] = i
            if fields is None:
                self.fields[i, :] = rotations.vec_random() * field_disorder
            else:
                self.fields[i, :] = fields[i, :]


    def getEnergy(self, coords):
        """
        coords is a list of (theta, phi) spherical coordinates of the spins
        where phi is the azimuthal angle (angle to the z axis) 
        """
        coords3 = coords2ToCoords3(coords)

        E = 0.
        for edge in self.G.edges():
            u = self.indices[edge[0]]
            v = self.indices[edge[1]]
            E -= np.dot(coords3[u, :], coords3[v, :])

        Efields = -np.sum(self.fields * coords3)

        return E + Efields

    def getEnergyGradient(self, coords):
        """
        coords is a list of (theta, phi) spherical coordinates of the spins
        where phi is the azimuthal angle (angle to the z axis) 
        """
        coords3 = coords2ToCoords3(coords)
        coords2 = coords

        E = 0.
        grad3 = np.zeros([self.nspins, 3])
        for edge in self.G.edges():
            u = self.indices[edge[0]]
            v = self.indices[edge[1]]
            E -= np.dot(coords3[u, :], coords3[v, :])

            grad3[u, :] -= coords3[v, :]
            grad3[v, :] -= coords3[u, :]

        Efields = -np.sum(self.fields * coords3)
        grad3 -= self.fields

        grad2 = grad3ToGrad2(coords2, grad3)
        grad2 = np.reshape(grad2, self.nspins * 2)

        return E + Efields, grad2


def normalize_spins(v3):
    v = v3.reshape([-1, 3])
    norms = np.sqrt((v * v).sum(1))
    v = v / norms[:, np.newaxis]
    v = v.reshape(-1)
    v3[:] = v[:]


def test_basin_hopping(pot, angles):  # pragma: no cover
    from pele.basinhopping import BasinHopping
    from pele.takestep.displace import RandomDisplacement
    from pele.takestep.adaptive import AdaptiveStepsize

    takestep = RandomDisplacement(stepsize=np.pi / 4)
    takestepa = AdaptiveStepsize(takestep, frequency=20)

    bh = BasinHopping(angles, pot, takestepa, temperature=1.01)
    bh.run(20)

# def test():
#    pi = np.pi
#    L = 8
#    nspins = L**2
#    
#    #phases = np.zeros(nspins)
#    pot = HeisenbergModel( dim = [L,L], field_disorder = 1.) #, phases=phases)
#    
#    coords = np.zeros([nspins, 2])
#    for i in range(nspins): 
#        vec = rotations.vec_random()
#        coords[i,:] = make2dVector(vec)
#    coords = np.reshape(coords, [nspins*2])
#    if False:
#        normfields = np.copy(pot.fields)
#        for i in range(nspins): normfields[i,:] /= np.linalg.norm(normfields[i,:])
#        coords = coords3ToCoords2( np.reshape(normfields, [nspins*3] ) )
#        coords  = np.reshape(coords, nspins*2)
#    #print np.shape(coords)
#    coordsinit = np.copy(coords)
#    
#    #print "fields", pot.fields
#    print coords
#    
#    if False:
#        coords3 = coords2ToCoords3(coords)
#        coords2 = coords3ToCoords2(coords3)
#        print np.reshape(coords, [nspins,2])
#        print coords2
#        coords3new = coords2ToCoords3(coords2)
#        print coords3
#        print coords3new
#
#    e = pot.getEnergy(coords)
#    print "energy ", e
#    if False:
#        print "numerical gradient"
#        ret = pot.getEnergyGradientNumerical(coords)
#        print ret[1]
#        if True:
#            print "analytical gradient"
#            ret2 = pot.getEnergyGradient(coords)
#            print ret2[1]
#            print ret[0]
#            print ret2[0]
#            print "ratio"
#            print ret2[1] / ret[1]
#            print "inverse sin"
#            print 1./sin(coords)
#            print cos(coords)
#
#    
#    print "try a quench"
#    from pele.optimize import mylbfgs
#    ret = mylbfgs(coords, pot, iprint=1)
#    
#    print "quenched e = ", ret.energy, "funcalls", ret.nfev
#    print ret.coords
#    with open("out.spins", "w") as fout:
#        s = coords2ToCoords3( ret.coords )
#        h = pot.fields
#        c = coords2ToCoords3( coordsinit )
#        for node in pot.G.nodes():
#            i = pot.indices[node]
#            fout.write( "%g %g %g %g %g %g %g %g %g %g %g\n" % (node[0], node[1], \
#                s[i,0], s[i,1], s[i,2], h[i,0], h[i,1], h[i,2], c[i,0], c[i,1], c[i,2] ) )
#    
#    coords3 = coords2ToCoords3( ret.coords )
#    m = np.linalg.norm( coords3.sum(0) ) / nspins
#    print "magnetization after quench", m
#    
#    test_basin_hopping(pot, coords)
#
#def test_potential():
#    L=4
#    nspins=L*L
#    pot = HeisenbergModel( dim = [L,L], field_disorder = 1.) #, phases=phases)
#    
#    coords = np.zeros([nspins, 2])
#    for i in range(nspins): 
#        vec = rotations.vec_random()
#        coords[i,:] = make2dVector(vec)
#    coords = np.reshape(coords, [nspins*2])
#    
#    coords[1] = 2.*np.pi - coords[1]
#    
#    pot.test_potential(coords)
#
#def test_potential_constrained():
#    L=4
#    nspins=L*L
#    pot = HeisenbergModelConstraint( dim = [L,L], field_disorder = 1.) #, phases=phases)
#    
#    coords = np.zeros([nspins, 3])
#    for i in range(nspins): 
#        coords[i,:] = rotations.vec_random()
#    coords = coords.reshape(-1)
#    
#    pot.test_potential(coords)
#
#
#if __name__ == "__main__":
#    test_potential_constrained()
##    test()

