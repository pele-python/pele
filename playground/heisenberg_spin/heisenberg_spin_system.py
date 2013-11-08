import numpy as np

from pele.potentials import HeisenbergModelRA
from pele.potentials import heisenberg_spin as hs
from pele.systems import BaseSystem
from pele.utils.rotations import vec_random
from pele.systems._opengl_tools import draw_cone
from pele.utils import rotations

def align_vectors_rmat(v1, v2):
    """return the rotation matrix that will align vector 2 with vector 1"""
    vx = np.cross(v1, v2)
    theta = np.arccos(np.dot(v1, v2))
    aa = theta * vx / np.linalg.norm(vx)
    return rotations.aa2mx(aa)


class spin3d_mindist(object):
    def __init__(self, pot):
        self.pot = pot
    
    def __call__(self, xa, xb):
        sa = hs.coords2ToCoords3(xa).reshape([-1,3])
        sb = hs.coords2ToCoords3(xb).reshape([-1,3])
        meansa = sa.sum(0)
        meansb = sb.sum(0)
        
        meansa /= np.linalg.norm(meansa)
        meansa /= np.linalg.norm(meansb)
        
        rmat = align_vectors_rmat(meansa, meansb)
        
        for i in xrange(self.pot.npins):
            sb[i,:] = np.dot(rmat, sb[i,:])
        
        dist = np.linalg.norm(sb.reshape(-1) - sa.reshape(-1))
        
        xb = hs.coords3ToCoords2(sb)
        return dist, xa, xb        
        
def spin3d_mindist_norot(xa, xb):
    sa = hs.coords2ToCoords3(xa)
    sb = hs.coords2ToCoords3(xb)
    dist = np.linalg.norm(sb - sa)
    return dist, xa, xb

        

class HeisenbergSystem(BaseSystem):
    def __init__(self, dims=[4,4], field_disorder=1.):
        BaseSystem.__init__(self)
        self.dims = dims
        self.field_disorder = field_disorder
        self.nspins = np.prod(dims)
        self.pot = self.get_potential()
        self.params.takestep.stepsize = np.pi / 4
        self.params.takestep.verbose = True
        
    def node2xyz(self, node):
        return np.array([float(x) for x in [node[0], node[1], 0]])
    
    def get_potential(self):
        try:
            return self.pot
        except AttributeError:
            return HeisenbergModelRA(dim=self.dims, field_disorder=self.field_disorder)
    
    def get_random_configuration(self):
        coords = np.zeros([self.nspins, 2])
        for i in range(self.nspins):
            vec = vec_random()
            coords[i,:] = hs.make2dVector(vec)
        return coords.reshape(-1)
    
    def get_mindist(self):
        # minimize the overall angle between the spins
        return spin3d_mindist_norot
    
    def get_orthogonalize_to_zero_eigenvectors(self):
        return None
    
    def get_metric_tensor(self, coords):
        return None

    def draw(self, coords, index):
        d = .4
        r = .04
        com = sum(self.node2xyz(node) for node in self.pot.G.nodes())
        com /= self.pot.nspins
        coords = coords.reshape([-1,2])
        for node in self.pot.G.nodes():
            xyz = self.node2xyz(node)
            spin2 = coords[self.pot.indices[node],:]
            spin3 = hs.make3dVector(spin2)
            x1 = xyz - 0.5 * spin3 * d - com
            x2 = xyz + 0.5 * spin3 * d - com
            draw_cone(x1, x2, rbase=r)
    
        

    
def rungui():
    from pele.gui import run_gui
    system = HeisenbergSystem(field_disorder=1.)
#    bh = system.get_basinhopping()
#    bh.run(10)
    run_gui(system)

def test_pot():
    np.random.seed(0)
    from pele.gui import run_gui
    system = HeisenbergSystem(field_disorder=4.)
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
#    test_pot()
    rungui()