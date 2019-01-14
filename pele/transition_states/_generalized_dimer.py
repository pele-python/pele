from __future__ import print_function
import numpy as np

from pele.transition_states import FindLowestEigenVector
from pele.optimize import Result
from pele.utils import rotations
from pele.transition_states._generalized_hef import _HybridEigenvectorWalker
from pele.transition_states._dimer_translator import _DimerTranslator


class GeneralizedDimer(object):
    """Use the generalized dimer method to find a saddle point
    
    This should be considered experimental.  It works, but I haven't spent
    enough time on it to make it very robust and efficient.  I would recommend
    using FindTransitionState instead.
    
    Parameters
    ----------
    coords : array
        The starting point for the optimization
    potential : the potential object
    eigenvec0 : array
        An initial guess for the smallest eigenvector (dimer direction)
    rotational_steps : int
        The number of iterations for optimizing the smallest eigenvector
        (rotating the dimer) between translational moves
    translational_steps : int
        The number of translational steps before re-optimizing the 
        smallest eigenvector 
    maxiter : int
        the maximum number of iterations
    minimizer_kwargs : dict
        these keyword arguments are passed to the minimizer that does the 
        translational steps
    leig_kwargs : dict
        these keyword arguments are passed to the class for optimizing 
        the smallest eigenvector
        
    Notes
    -----
    The dimer method defines a way of walking along an energy surface
    towards a saddle point.  The dimer is defined by a location and 
    a direction, vec.  vec will point along the direction of largest
    negative curvature, which corresponds to the eigenvector of the Hessian
    matrix with the smallest eigenvalue.  The eigenvalue is the 
    curvature.  The eigenvector can be computed analytically if the 
    Hessian is available, but this is often computationally expensive
    normally you would use gradients and an optimization function to 
    converge towards the direction of largest curvature.  If the
    gradient is inverted along the direction of the largest negative
    curvature then following that gradient will lead you towards a
    saddle point.  As you walk towards the saddle point the direction
    of largest negative curvature will change, so it will need to be 
    periodically re-optimized.  Thus the main iteration loop is
    
    1. Optimize the rotation of the dimer so that it is aligned with
    the direction of largest negative curvature.
    
    2. Translate the dimer uphill following the direction of largest
    negative curvature while simultaneously walking downhill in all
    perpendicular directions.
    
    3. If converged, then end, else go to 1.   
    
    """

    def __init__(self, coords, potential, eigenvec0=None,
                 rotational_steps=20,
                 translational_steps=10,
                 maxiter=500,
                 leig_kwargs=None,
                 translator_kwargs=None,
                 dimer=True,
    ):
        coords = coords.copy()
        self.rotational_steps = rotational_steps
        self.translational_steps = translational_steps
        self.maxiter = maxiter
        self.iter_number = 0

        # check the keyword dictionaries
        if translator_kwargs is None: translator_kwargs = {}
        if leig_kwargs is None: leig_kwargs = {}

        # set up the initial guess for the eigenvector
        if eigenvec0 is None:
            eigenvec0 = rotations.vec_random_ndim(coords.shape)
        eigenvec0 /= np.linalg.norm(eigenvec0)
        assert coords.shape == eigenvec0.shape

        # set up the object that will maintain the rotation of the dimer
        self.rotator = FindLowestEigenVector(coords, potential, eigenvec0=eigenvec0, **leig_kwargs)

        # set up the object that will translate the dimer
        if dimer:
            self.translator = _DimerTranslator(coords, potential, eigenvec0, **translator_kwargs)
        else:
            self.translator = _HybridEigenvectorWalker(coords, potential, eigenvec0, **translator_kwargs)

    def get_true_energy_gradient(self, coords):
        """return the true energy and gradient"""
        return self.translator.get_true_energy_gradient(coords)


    # def get_true_energy(self):
    # """return the true energy"""
    # return self.translator.get_energy()
    #
    # def get_true_gradient(self):
    # """return the true gradient"""
    # # these are stored in dimer_potential
    # return self.translator.get_gradient()

    def get_coords(self):
        """return the current location of the dimer"""
        mret = self.translator.get_result()
        return mret.coords

    def stop_criterion_satisfied(self):
        """return True if the stop criterion is satisfied"""
        return self.translator.stop_criterion_satisfied() and self.rotator.stop_criterion_satisfied()

    def one_iteration(self):
        """do one iteration"""
        # update the eigenvector (rotate the dimer)
        # print "rotating dimer"
        energy, gradient = self.get_true_energy_gradient(self.get_coords())
        self.rotator.update_coords(self.get_coords(), gradient=gradient)
        ret = self.rotator.run(self.rotational_steps)

        # update the eigenvector and eigenvalue in the dimer_potential
        self.translator.update_eigenvec(ret.eigenvec, ret.eigenval)

        # translate the dimer
        # print "translating dimer"
        self.translator.run(self.translational_steps)

        self.iter_number += 1


    def run(self):
        """the main iteration loop"""
        while not self.stop_criterion_satisfied() and self.iter_number < self.maxiter:
            self.one_iteration()

        return self.get_result()


    def get_result(self):
        """return a results object"""
        trans_res = self.translator.get_result()

        rot_result = self.rotator.get_result()

        res = Result()
        res.eigenval = rot_result.eigenval
        res.eigenvec = rot_result.eigenvec
        res.coords = trans_res.coords
        res.energy, res.grad = self.get_true_energy_gradient(res.coords)
        res.rms = trans_res.rms
        res.nfev = rot_result.nfev + trans_res.nfev
        res.nsteps = self.iter_number
        res.success = self.stop_criterion_satisfied()

        if res.eigenval > 0:
            res.success = False

        return res


#
# testing only below here
#

class PotWrapper(object):  # pragma: no cover
    def __init__(self, pot):
        self.pot = pot
        self.nfev = 0

    def getEnergyGradient(self, x):
        self.nfev += 1
        return self.pot.getEnergyGradient(x)


def compare_HEF(x0, evec0, system, **kwargs):  # pragma: no cover
    from pele.transition_states import findTransitionState

    pot = PotWrapper(system.get_potential())
    ret = findTransitionState(x0, pot, eigenvec0=evec0, orthogZeroEigs=None, **kwargs)
    print(ret.eigenval)
    print(pot.nfev)
    print(ret.rms)


def get_x0():  # pragma: no cover
    from pele.systems import LJCluster

    natoms = 31
    system = LJCluster(natoms)
    db = system.create_database()
    bh = system.get_basinhopping(db, outstream=None)
    while db.number_of_minima() < 2:
        bh.run(1)

    mindist = system.get_mindist()
    m1, m2 = db.minima()[:4]
    d, x1, x2 = mindist(m1.coords, m2.coords)

    x0 = (x1 + x2) / 2
    evec0 = x2 - x1

    return system, x0, evec0


def test():  # pragma: no cover
    system, x0, evec0 = get_x0()

    dimer = GeneralizedDimer(x0.copy(), system.get_potential(),
                             eigenvec0=evec0,
                             # dimer_kwargs=dict(n_translational_steps=5,
                             # n_rotational_steps=20),
                             # minimizer_kwargs=dict(iprint=1, tol=4e-5, nsteps=2000)
    )
    ret = dimer.run()

    print("eigenvalue", ret.eigenval)
    print("function evaluations", ret.nfev)
    print("rms", ret.rms)

    print("\n\nnow with hybrid eigenvector following")
    compare_HEF(x0, evec0, system, nsteps_tangent1=5, nsteps_tangent2=5,
                lowestEigenvectorQuenchParams={"nsteps": 20, "tol": 1e-4}, iprint=10)


if __name__ == "__main__":
    test()

