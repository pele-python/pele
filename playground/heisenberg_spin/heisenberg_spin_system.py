import numpy as np

from pele.potentials import HeisenbergModelRA, HeisenbergModel
from pele.potentials import heisenberg_spin as hs
from pele.systems import BaseSystem
from pele.utils.rotations import vec_random
from pele.systems._opengl_tools import draw_cone
from pele.utils import rotations
from pele.landscape import smooth_path
from pele.utils.frozen_atoms import FrozenCoordsConverter, FrozenPotWrapper
#from pele.transition_states import InterpolatedPath

#def align_vectors_rmat(v1, v2):
#    """return the rotation matrix that will align vector 2 with vector 1"""
#    vx = np.cross(v1, v2)
#    theta = np.arccos(np.dot(v1, v2))
#    aa = theta * vx / np.linalg.norm(vx)
#    return rotations.aa2mx(aa)

def interpolate_spin(v1, v2, t):
    vx = np.cross(v2, v1)
    theta = np.arccos(np.dot(v1, v2))
    theta *= (1.0 - t)
    aa = theta * vx / np.linalg.norm(vx)
    mx = rotations.aa2mx(aa)
    v3 = np.dot(mx, v2)
    return v3 / np.linalg.norm(v3)

def interpolate_spins(initial, final, t, i3=None, f3=None):
    if i3 is None:
        i3 = hs.coords2ToCoords3(initial)
    if f3 is None:
        f3 = hs.coords2ToCoords3(final)
    nspins = i3.size / 3
    i3 = i3.reshape([-1,3])
    f3 = f3.reshape([-1,3])
    xnew = np.zeros(i3.shape)
    for i in xrange(nspins):
        xnew[i,:] = interpolate_spin(i3[i,:], f3[i,:], t)
    return hs.coords3ToCoords2(xnew.reshape(-1)).reshape(-1)

#def smooth_path_spins(path, ):
#    return smooth_path(path, mindist, density, interpolator)
#    i3 = hs.coords2ToCoords3(initial)
#    f3 = hs.coords2ToCoords3(final)
    

#class spin3d_mindist(object):
#    def __init__(self, pot):
#        self.pot = pot
#    
#    def __call__(self, xa, xb):
#        sa = hs.coords2ToCoords3(xa).reshape([-1,3])
#        sb = hs.coords2ToCoords3(xb).reshape([-1,3])
#        meansa = sa.sum(0)
#        meansb = sb.sum(0)
#        
#        meansa /= np.linalg.norm(meansa)
#        meansa /= np.linalg.norm(meansb)
#        
#        rmat = align_vectors_rmat(meansa, meansb)
#        
#        for i in xrange(self.pot.npins):
#            sb[i,:] = np.dot(rmat, sb[i,:])
#        
#        dist = np.linalg.norm(sb.reshape(-1) - sa.reshape(-1))
#        
#        xb = hs.coords3ToCoords2(sb)
#        return dist, xa, xb        


def spin3d_mindist_norot(xa, xb):
    """
    would be better to use the spherical law of cosines to compute
    the angle between the vectors
    http://en.wikipedia.org/wiki/Great-circle_distance
    """
    sa = hs.coords2ToCoords3(xa).reshape([-1,3])
    sb = hs.coords2ToCoords3(xb).reshape([-1,3])
    # get an array of the dot product between the spins
    dots = np.sum(sa * sb, axis=1)
    angles = np.arccos(dots)
    dist = angles.sum() / np.pi
#    dist = np.linalg.norm(sb - sa)
    return dist, xa, xb

def normalize_2dspins(coords):
    """there is probably a better way to do this."""
    return hs.make2dVector(hs.make3dVector(coords))

def spin3d_distance(xa, xb, distance=True, grad=True):
    """
    the coorindates are in the form (theta, phi) = xa[0,:]
    where phi is the polar angle (the angle to the z axis)
    and theta is the azimulth angle to the x-axis (I think... it could be y axis).
    1) compute delta phi and delta theta
    2) apply symmetries until -pi/2 < delta phi < pi/2
    3) apply symmetries until -pi < delta theta < pi
    4) normalize delta theta by sin(phi) 
    """
    import sys
    sys.stderr.write("Warning, this distance function (spin3d_distance) is wrong.  The NEB will be messed up\n")
    from pele.angleaxis import _aadist
    dist, temp1, temp2 = spin3d_mindist_norot(xa, xb)
    S = np.eye(3)
    print "THIS IS WHAT I WAS DOING BEFORE I WENT TO THE PUB"
    xa = normalize_2dspins(xa)
    xb = normalize_2dspins(xb)
    grad = xa - xb
    if distance:
        dist, xa, xb = spin3d_mindist_norot(xa, xb)
    return dist, grad

class HeisenbergSystem(BaseSystem):
    def __init__(self, dims=[4,4], field_disorder=1., disorder=False):
        BaseSystem.__init__(self)
        self.dims = dims
        self.field_disorder = field_disorder
        self.nspins = np.prod(dims)

        self.one_frozen = False
        if field_disorder == 0. or not disorder:
            self.one_frozen = True

        self.pot = self.get_potential()
        
        self.setup_params(self.params)
    
    def setup_params(self, params):
        params.takestep.stepsize = np.pi #/ 2.
        params.takestep.verbose = True
        params.double_ended_connect.local_connect_params.NEBparams.interpolator = interpolate_spins
        params.double_ended_connect.local_connect_params.NEBparams.image_density = 10
        params.double_ended_connect.local_connect_params.NEBparams.reinterpolate = 50
        params.double_ended_connect.local_connect_params.NEBparams.distance = spin3d_distance
        params.structural_quench_params.tol = 1e-6
                
    def node2xyz(self, node):
        return np.array([float(x) for x in [node[0], node[1], 0]])
    
    def get_potential(self):
        try:
            return self.pot
        except AttributeError:
            if self.one_frozen:
                # freeze one spin to remove the global rotational symmetry
                base_pot = HeisenbergModel(dim=self.dims, field_disorder=0.)
                reference_coords = np.zeros(base_pot.nspins * 2)
                n = reference_coords.size
                frozen_dof = np.array([n-2, n-1])
#                self.coords_converter = FrozenCoordsConverter(reference_coords, frozen_dof)
                self.pot = FrozenPotWrapper(base_pot, reference_coords, frozen_dof)
                self.coords_converter = self.pot.coords_converter
                self.pot.G = base_pot.G
                self.pot.indices = base_pot.indices
                return self.pot
            else:
                return HeisenbergModelRA(dim=self.dims, field_disorder=self.field_disorder)
    
    def get_random_configuration(self):
        nspins = self.nspins
        if self.one_frozen:
            nspins -= 1
        coords = np.zeros([nspins, 2])
        for i in range(nspins):
            vec = vec_random()
            coords[i,:] = hs.make2dVector(vec)
        return coords.reshape(-1)
    
    def get_mindist(self):
        # minimize the overall angle between the spins
        return spin3d_mindist_norot
    
    def get_orthogonalize_to_zero_eigenvectors(self):
        return None
    
    def get_metric_tensor(self, coords):
        """
        metric tensor on surface of sphere
        
        g = [[sin(phi)**2  0
              0              1]]
        
        where phi is angle measured from the z-axis.  Which is coords[:,1].  
        The distance vector is then written
        
        ds**2 = sin(phi)**2 * dtheta**2 + dphi**2
        """
#        return None
        coords = coords.reshape([-1,2])
        sinphi2 = np.sin(coords[:,1])**2
        n = coords.size
        g = np.eye(n).reshape(-1)
        # set every second diagonal component to sintheta
#        g[n+1::2*(n+1)] = sintheta2
        g[::2*(n+1)] = sinphi2
        return g.reshape([n,n])

    def draw(self, coords, index):
        if self.one_frozen:
            coords = self.coords_converter.get_full_coords(coords)
        d = .4
        r = .04
        nspins = coords.size / 2
        com = sum(self.node2xyz(node) for node in self.pot.G.nodes())
        com /= nspins
        coords = coords.reshape([-1,2])
        for node in self.pot.G.nodes():
            xyz = self.node2xyz(node)
            spin2 = coords[self.pot.indices[node],:]
            spin3 = hs.make3dVector(spin2)
            x1 = xyz - 0.5 * spin3 * d - com
            x2 = xyz + 0.5 * spin3 * d - com
            draw_cone(x1, x2, rbase=r)
    
    def smooth_path(self, images, **kwargs):
        mindist = self.get_mindist()
        path = smooth_path(images, mindist, interpolator=interpolate_spins)
        return path
        

    
def rungui():
    from pele.gui import run_gui
    system = HeisenbergSystem(field_disorder=3.1, disorder=True)
#    system = HeisenbergSystem(field_disorder=0., disorder=False)

#    x = system.get_random_configuration()
#    print system.get_metric_tensor(x)
#    return
#    bh = system.get_basinhopping()
#    bh.run(10)
    run_gui(system)

def test_eigs():
    from pele.transition_states import findLowestEigenVector
    from pele.thermodynamics import normalmodes
    system = HeisenbergSystem(field_disorder=3.1, disorder=True)
    pot = system.get_potential()
    x = system.get_random_configuration()
    x = system.get_random_minimized_configuration().coords
    ret = findLowestEigenVector(x, pot, orthogZeroEigs=None)
    hess = pot.getHessian(x)
    freq, modes = normalmodes(hess, metric=None)
    print "lowest eig from hess", freq[0]
    print "lowest eig from RR  ", ret.eigenval
    

def test_pot():
    np.random.seed(0)
    from pele.gui import run_gui
    system = HeisenbergSystem(field_disorder=3., disorder=True)
    pot = system.get_potential()
    x = system.get_random_configuration()
    opt = system.get_minimizer(iprint=1)
    print opt
    from pele.optimize import lbfgs_py as optimizer
    from pele.optimize import lbfgs_cpp as optimizer
    from pele.optimize import mylbfgs as optimizer
    opt = lambda coords: optimizer(coords, system.get_potential(), iprint=1)
    ret = opt(x)
    xnew = ret.coords
    print ret.energy
    print ret.grad
    print ret
#    print ret.coords
    print "computed grad", pot.getGradient(xnew)
    
    
    pot.test_potential(xnew)
    
        
if __name__ == "__main__":
    print np.arccos(-.99)
    print np.arccos(0.99)
    print np.arccos(0.0)
    np.random.seed(0)
#    test_pot()
#    test_eigs()
    rungui()
