import numpy as np

from pele.potentials import XYModel
from pele.potentials.xyspin import angle_to_2dvector
from pele.systems import BaseSystem
from pele.landscape import smooth_path
from pele.utils.frozen_atoms import FrozenPotWrapper

def normalize_spins(x):
    L = 2. * np.pi
    x -= L * np.floor(x / L)
    return x

def spin_distance_1d(x1, x2):
    dx = x1 - x2
    # apply periodic boundary conditions
    L = 2. * np.pi
    dx -=  L * np.round(dx / L)
    return np.linalg.norm(dx)

def spin_mindist_1d(x1, x2):
    # apply periodic boundary conditions
    L = 2. * np.pi
    offset = L * np.round((x1 - x2) / L)
    x2 += offset
    assert np.max(np.abs(x1-x2)) <= L/2.
    return np.linalg.norm(x1-x2), x1, x2

class XYModlelSystem(BaseSystem):
    def __init__(self, dim=[4, 4], phi_disorder=np.pi, phases=None):
        BaseSystem.__init__(self)
        self.one_frozen = True
        self.dim = dim
        self.phi_disorder = phi_disorder
        self.pot = self.get_potential(phases=phases)
        self.nspins = np.prod(dim)
        
        self.setup_params(self.params)

    def setup_params(self, params):
        params.takestep.stepsize = np.pi# / 2.
        params.takestep.verbose = True
#        params.double_ended_connect.local_connect_params.NEBparams.interpolator = interpolate_spins
        params.double_ended_connect.local_connect_params.NEBparams.image_density = .8
        params.double_ended_connect.local_connect_params.NEBparams.iter_density = 50.
        params.double_ended_connect.local_connect_params.NEBparams.reinterpolate = 50
        params.double_ended_connect.local_connect_params.NEBparams.adaptive_nimages = True
        params.double_ended_connect.local_connect_params.NEBparams.adaptive_niter = False
#        params.double_ended_connect.local_connect_params.NEBparams.distance = spin3d_distance
        params.structural_quench_params.tol = 1e-6
        params.database.overwrite_properties = False
        
        params.basinhopping.insert_rejected = True

    def get_system_properties(self):
        return dict(potential="XY spin glass",
                    phases=self.pot.get_phases(),
                    dim=self.dim,
                    Lx=self.dim[0],
                    Ly=self.dim[1],
                    )

    def get_potential(self, phases=None):
        try:
            return self.pot
        except AttributeError:
            assert self.one_frozen
            base_pot = XYModel(dim=self.dim, phi=self.phi_disorder, phases=phases)
            reference_coords = np.zeros(base_pot.nspins)
            n = reference_coords.size
            frozen_node = (0,0)
            frozen_index = base_pot.indices[frozen_node]
            frozen_dof = np.array([frozen_index])
            print "making frozen spin at", self.node2xyz(base_pot.index2node[frozen_index])
#                self.coords_converter = FrozenCoordsConverter(reference_coords, frozen_dof)
            self.pot = FrozenPotWrapper(base_pot, reference_coords, frozen_dof)
#            self.coords_converter = self.pot.coords_converter
#            self.pot.G = base_pot.G
#            self.pot.indices = base_pot.indices
            return self.pot

    def get_orthogonalize_to_zero_eigenvectors(self):
        return None
    
    def get_metric_tensor(self, coords):
        return None
    
    def get_nzero_modes(self):
        if self.one_frozen:
            return 0
        else:
            return 1

    def get_pgorder(self, coords):
        return 1
    
    def get_mindist(self):
        assert self.one_frozen
        return spin_mindist_1d

    def get_compare_exact(self):
        mindist = self.get_mindist()
        return lambda x1, x2: mindist(x1, x2)[0] < .1

    def smooth_path(self, path, **kwargs):
        mindist = self.get_mindist()
        return smooth_path(path, mindist, **kwargs)

    def get_random_configuration(self):
        if self.one_frozen:
            n = self.nspins - 1
        else:
            n = self.nspins
        return np.random.uniform(0, 2.*np.pi, n)

    def node2xyz(self, node):
        return np.array([float(x) for x in [node[0], node[1], 0]])
    
    def create_database(self, *args, **kwargs):
        return BaseSystem.create_database(self, *args, **kwargs)

    def get_takestep(self, **kwargs):
        """return the takestep object for use in basinhopping, etc.
        
        default is random displacement with adaptive step size 
        adaptive temperature
        
        See Also
        --------
        pele.takestep
        """
        if self.phi_disorder > 0.01:
            return super(XYModlelSystem, self).get_takestep(**kwargs)
        # if no disorder, turn off adaptive step and temperature.
        from pele.takestep import RandomDisplacement
        kwargs = dict(self.params["takestep"].items() + kwargs.items())
        try:
            stepsize = kwargs.pop("stepsize")
        except KeyError:
            stepsize = 2.*np.pi
        takeStep = RandomDisplacement(stepsize=stepsize)
        return takeStep

    def draw(self, coords, index):
        from pele.systems._opengl_tools import draw_cone
#        if self.one_frozen:
#            coords = self.coords_converter.get_full_coords(coords)
        d = .4
        r = .04
        if self.one_frozen:
            coords = self.pot.coords_converter.get_full_coords(coords)
        nspins = coords.size
        com = sum(self.node2xyz(node) for node in self.pot.G.nodes())
        com /= nspins
        for node in self.pot.G.nodes():
            xyz = self.node2xyz(node)
            i = self.pot.indices[node]
            spin2 = angle_to_2dvector(coords[i])
            spin3 = np.zeros(3)
            spin3[:2] = spin2
            x1 = xyz - 0.5 * spin3 * d - com
            x2 = xyz + 0.5 * spin3 * d - com
            draw_cone(x1, x2, rbase=r)
    
def normalize_spins_db(db):
    for m in db.minima():
        x = normalize_spins(m.coords)
        print np.max(x), np.min(x)
        m.coords = x
    db.session.commit()
    
def run_gui():
    from pele.gui import run_gui
    system = XYModlelSystem(dim=[24,24], phi_disorder=np.pi)
    run_gui(system)

def run_gui_db(dbname="xy_10x10.sqlite"):
    from pele.gui import run_gui
    from pele.storage import Database
    try:
        db = Database(dbname, createdb=False)
        phases = db.get_property("phases").value()
    except IOError:
        phases=None
    system = XYModlelSystem(dim=[10,10], phi_disorder=np.pi, phases=phases)
    run_gui(system, db=dbname)

def run_gui_nodisorder(L=24):
    from pele.gui import run_gui
    dim=[L,L]
    system = XYModlelSystem(dim=dim, phi_disorder=0.)
    system.params.basinhopping.temperature=10.
    dbname="xy_%dx%d_nodisorder.sqlite" %(L,L)
    run_gui(system, db=dbname)

def test_potential():
    system = XYModlelSystem(dim=[10,10], phi_disorder=np.pi)
    pot = system.get_potential()
    pot.test_potential(system.get_random_configuration())


if __name__ == "__main__":
#    test_potential()
#    from pele.storage import Database
#    db = Database("20x20_no_disorder.sqlite")
#    normalize_spins_db(db)
#    run_gui_db()
    run_gui_nodisorder()
        
    
    
