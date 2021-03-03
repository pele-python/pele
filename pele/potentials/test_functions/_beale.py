from __future__ import print_function
import numpy as np

from pele.potentials import BasePotential
from pele.systems import BaseSystem


def makeplot2d(f, nx=100, xmin=None, xmax=None, zlim=None, show=True):  # pragma: no cover
    from matplotlib import cm
    import matplotlib.pyplot as plt
    import numpy as np

    ny = nx
    if xmin is None:
        xmin = f.xmin[:2]
    xmin, ymin = xmin
    if xmax is None:
        xmax = f.xmax[:2]
    xmax, ymax = xmax
    x = np.arange(xmin, xmax, (xmax - xmin) / nx)
    y = np.arange(ymin, ymax, (ymax - ymin) / ny)
    X, Y = np.meshgrid(x, y)
    Z = np.zeros(X.shape)
    for i in range(x.size):
        for j in range(y.size):
            xy = np.array([X[i, j], Y[i, j]])
            Z[i, j] = f.getEnergy(xy)

    if zlim is not None:
        Z = np.where(Z > zlim[1], -1, Z)

    fig = plt.figure()
    # ax = fig.gca(projection='3d')
    mesh = plt.pcolormesh(X, Y, Z, cmap=cm.coolwarm)
    # surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
    # linewidth=0, antialiased=False)

    # if zlim is not None:
    # ax.set_zlim(zlim)

    # ax.zaxis.set_major_locator(LinearLocator(10))
    # ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    fig.colorbar(mesh)  # surf, shrink=0.5, aspect=5)

    if show:
        plt.show()
    return fig


class Beale(BasePotential):
    target_E = 0.
    target_coords = np.array([3., 0.5])
    xmin = np.array([-4.5, -4.5])
    # xmin = np.array([0., 0.])
    xmax = np.array([4.5, 4.5])

    def getEnergy(self, coords):
        x, y = coords
        return (1.5 - x + x * y) ** 2 + (2.25 - x + x * y ** 2) ** 2 + (2.625 - x + x * y ** 3) ** 2

    def getEnergyGradient(self, coords):
        E = self.getEnergy(coords)
        x, y = coords
        dx = 2. * (1.5 - x + x*y) * (-1. + y) + 2. * (2.25 - x + x * y**2) * (-1. + y**2) + 2. * (2.625 - x + x * y**3) * (-1. + y**3)
        dy = 2. * (1.5 - x + x*y) * x +       2. * (2.25 - x + x * y**2) * (2. * y * x) + 2. * (2.625 - x + x * y**3) * (3. * x * y**2)
        return E, np.array([dx, dy])


class BealeSystem(BaseSystem):
    def get_potential(self):
        return Beale()

    def get_random_configuration(self, eps=1e-3):
        pot = self.get_potential()
        xmin, xmax = pot.xmin, pot.xmax
        x = np.random.uniform(xmin[0] + eps, xmax[0] - eps)
        y = np.random.uniform(xmin[1] + eps, xmax[1] - eps)
        return np.array([x, y])


def add_minimizer(pot, ax, minimizer, x0, **kwargs):  # pragma: no cover
    xcb = []

    def cb(coords=None, energy=None, rms=None, **kwargs):
        xcb.append(coords.copy())
        print("energy", energy, rms)

    minimizer(x0, pot, events=[cb], **kwargs)
    xcb = np.array(xcb)
    ax.plot(xcb[:, 0], xcb[:, 1], '-o', label=minimizer.__name__)


def test_minimize():  # pragma: no cover
    # from pele.potentials.test_functions import BealeSystem
    # from pele.potentials.test_functions._beale import makeplot2d
    from pele.optimize import lbfgs_py, steepest_descent
    import matplotlib.pyplot as plt

    system = BealeSystem()
    pot = system.get_potential()
    fig = makeplot2d(pot, nx=60, show=False, zlim=[0, 50])
    ax = fig.gca()
    x0 = system.get_random_configuration(eps=.5)

    # add_minimizer(pot, ax, fire, x0, nsteps=200, maxstep=.1)
    add_minimizer(pot, ax, lbfgs_py, x0, nsteps=200, M=1)
    add_minimizer(pot, ax, steepest_descent, x0, nsteps=10000, dx=1e-4)
    plt.legend()


    # lbfgs_py(system.get_random_configuration(), pot, events=[callback], nsteps=100, M=1)
    plt.show()


def test1():  # pragma: no cover
    s = BealeSystem()
    f = s.get_potential()
    f.test_potential(f.target_coords)
    print("")
    f.test_potential(s.get_random_configuration())
    f.test_potential(np.array([1., 1.]))  # , print_grads=True)

    # from base_function import makeplot2d
    v = 3.
    makeplot2d(f, nx=60, zlim=[0, 100])


if __name__ == "__main__":
    test_minimize()

