import numpy as np
import os.path
import copy
import logging

from pele.optimize import Result
from pele.optimize import mylbfgs

__all__ = ["NEB", ]

logger = logging.getLogger("pele.connect.neb")

try:
    import scipy.linalg

    blas = lambda name, ndarray: scipy.linalg.get_blas_funcs((name,), (ndarray,))[0]
    bnorm = blas('nrm2', np.array([], dtype="float64"))

    def norm(x):
        assert x.dtype == "float64"
        return bnorm(x)
except Exception:
    norm = np.linalg.norm


def distance_cart(x1, x2, distance=True, grad=True):
    dist = None
    grad = x1 - x2
    if distance:
        dist = norm(grad) ** 2
    return dist, grad


class NEB(object):
    """Doubly nudged elastic band implementation

    Parameters
    ----------
    path: iteratable
        iteratable object of all images
    potential :
        the potential object
    distance : callable
        distance function for the elastic band
    k : float, optional
        the elastic band spring constant 
    adjustk_freq : integer
        frequency to adjust k, set to zero to disable
    adjustk_tol : float
        tolerance for adjusting k up or down
    adjustk_factor : float
        the multiplicative factor used to adjust k
    dneb: bool, optional
        do double nudging, default True
    with_springenergy: boolean, optional
        add the spring energy to the total energy of the band, default is False
    copy_potential : bool, optional
        if True a separate copy of the potential will be made for
        each image.  This can be used to keep neighbor lists from being rebuilt
        over and over again.
    quenchRoutine : callable
        the quench routine to use in optimizing the band
    quenchParams :
        parameters passed to the quench routine. The default tolerance for quench is 1e-4
    save_energies : bool
        if True and quenchParams.iprint exists and is positive, then the energy 
        along the NEB path will be printed to a file of the form "neb.EofS.####" 
        every iprint steps.
    verbose : integer
        verbosity level
    events : list of callables
        list of callback functions called just before getEnergyGradient returns
    use_minimizer_callback: boolean, optional
        use the callback function of theminimizer to adjust k, it is not recommended to change this to false

    Notes
    -----
    use the Interpolation tools in this package to construct the initial path
    from a starting and ending point
    
    See Also
    --------
    InterpolatedPath : used to construct required parameter `path`
    InterpolatedPathDensity : alternate interpolater
    """

    def __init__(self, path, potential, distance=distance_cart, k=100.0, adjustk_freq=0, adjustk_tol=0.1,
                 adjustk_factor=1.05, with_springenergy=False, dneb=True, copy_potential=False, quenchParams=None,
                 quenchRoutine=None, save_energies=False, verbose=-1, events=None, use_minimizer_callback=True):
        if quenchParams is None: quenchParams = dict()
        self.distance = distance
        self.potential = potential
        self.k = k
        self.verbose = verbose
        self.use_minimizer_callback = use_minimizer_callback

        if events is None:
            self.events = []
        else:
            self.events = events

        nimages = len(path)
        self.nimages = nimages
        assert self.nimages >= 3
        self.copy_potential = copy_potential
        if self.copy_potential:
            self.potential_list = [copy.deepcopy(self.potential) for _ in range(nimages)]

        self.getEnergyCount = 0
        self.printStateFile = None
        self.iprint = -1
        self.save_energies = save_energies

        self.quenchRoutine = quenchRoutine
        self.quenchParams = quenchParams.copy()

        if "tol" not in quenchParams:
            self.quenchParams["tol"] = 1e-4

        self.adjustk_freq = adjustk_freq
        self.adjustk_tol = adjustk_tol
        self.adjustk_factor = adjustk_factor

        # initialize coordinate&gradient array
        self.coords = np.zeros([nimages, path[0].size])
        self.energies = np.zeros(nimages)
        self.isclimbing = []
        for i in range(nimages):
            self.isclimbing.append(False)

        # copy initial path
        for i, x in zip(range(self.nimages), path):
            self.coords[i, :] = x

        for i in range(0, nimages):
            self.energies[i] = potential.getEnergy(self.coords[i, :])
        # the active range of the coords, endpoints are fixed
        self.active = self.coords[1:nimages - 1, :]

        self.dneb = dneb
        self.with_springenergy = with_springenergy

        self.distances = np.zeros(self.nimages - 1)
        self.rms = 0.

    def optimize(self, quenchRoutine=None,
                 **kwargs):
        """
        Optimize the band

        Notes
        -----
        the potential for the NEB optimization is not Hamiltonian.  This
        means that there is no meaningful energy associated with the
        potential.  Therefore, during the optimization, we can do gradient
        following, but we can't rely on the energy for, e.g. determining
        step size, which is the default behavior for many optimizers.  This
        can be worked around by choosing a small step size and a large
        maxErise, or by using an optimizer that uses only gradients.

        scipy.lbfgs_b seems to work with NEB pretty well, but lbfgs_py and
        mylbfgs tend to fail.  If you must use one of those try, e.g.
        maxErise = 1., maxstep=0.01, tol=1e-2

        :quenchRoutine: quench algorithm to use for optimization.

        :quenchParams: parameters for the quench """
        if quenchRoutine is None:
            if self.quenchRoutine is None:
                quenchRoutine = mylbfgs
            else:
                quenchRoutine = self.quenchRoutine
                # combine default and passed params.  passed params will overwrite default
        quenchParams = dict([("nsteps", 300)] +
                            list(self.quenchParams.items()) +
                            list(kwargs.items()))

        if "iprint" in quenchParams:
            self.iprint = quenchParams["iprint"]
        if "logger" not in quenchParams:
            quenchParams["logger"] = logging.getLogger("pele.connect.neb.quench")

        if self.use_minimizer_callback:
            quenchParams["events"] = [self._step]

        self.step = 0
        qres = quenchRoutine(
            self.active.reshape(self.active.size), self,
            **quenchParams)
        # if isinstance(qres, tuple): # for compatability with old and new quenchers
        # qres = qres[4]

        self.active[:, :] = qres.coords.reshape(self.active.shape)
        if self.copy_potential:
            for i in range(0, self.nimages):
                pot = self.potential_list[i]
                self.energies[i] = pot.getEnergy(self.coords[i, :])
        else:
            for i in range(0, self.nimages):
                self.energies[i] = self.potential.getEnergy(self.coords[i, :])

        res = Result()
        res.path = self.coords
        res.nsteps = qres.nsteps
        res.energy = self.energies
        res.rms = qres.rms
        res.success = False
        if qres.rms < quenchParams["tol"]:
            res.success = True

        return res

    def _getRealEnergyGradient(self, coordsall):
        # calculate real energy and gradient along the band. energy is needed for tangent
        # construction
        realgrad = np.zeros(coordsall.shape)
        if self.copy_potential:
            for i in range(1, self.nimages - 1):
                pot = self.potential_list[i]
                self.energies[i], realgrad[i, :] = pot.getEnergyGradient(coordsall[i, :])
        else:
            for i in range(1, self.nimages - 1):
                self.energies[i], realgrad[i, :] = self.potential.getEnergyGradient(coordsall[i, :])
        return realgrad

    def getEnergy(self, coords):
        """this is a very dumb way of getting the energy.  it should be replaced"""
        e, g = self.getEnergyGradient(coords)
        return e

    def getEnergyGradient(self, coords1d):
        """
        Calculates the gradient for the whole NEB. only use force based minimizer!

        coords1d:
            coordinates of the whole neb active images (no end points)
        """
        # make array access a bit simpler, create array which contains end images
        tmp = self.coords.copy()
        tmp[1:self.nimages - 1, :] = coords1d.reshape(self.active.shape)
        grad = np.zeros(self.active.shape)

        # calculate real energy and gradient along the band. energy is needed for tangent
        # construction
        realgrad = self._getRealEnergyGradient(tmp)

        # the total energy of images, band is neglected
        E = sum(self.energies)
        Eneb = 0
        # build forces for all images
        for i in range(1, self.nimages - 1):
            En, grad[i - 1, :] = self.NEBForce(
                self.isclimbing[i],
                [self.energies[i], tmp[i, :]],
                [self.energies[i - 1], tmp[i - 1, :]],
                [self.energies[i + 1], tmp[i + 1, :]],
                realgrad[i, :],
                i
            )
            Eneb += En
        if self.iprint > 0:
            if self.getEnergyCount % self.iprint == 0 and self.save_energies:
                self.printState()
        self.getEnergyCount += 1
        if not self.use_minimizer_callback:
            self._step(coords1d)

        rms = np.linalg.norm(grad) / np.sqrt(self.active.size)

        for event in self.events:
            event(path=tmp, energies=self.energies,
                  distances=self.distances, stepnum=self.getEnergyCount, rms=rms)

        return E + Eneb, grad.reshape(grad.size)

    def tangent_old(self, central, left, right, gleft, gright):
        """
        Old tangent construction based on average of neighbouring images

        coords1d:
            coordinates of the whole neb active images (no end points)
        """
        t = gleft / norm(gleft) - gright / norm(gright)
        return t / norm(t)

    def tangent(self, central, left, right, gleft, gright):
        """
        New uphill tangent formulation

        The method was  described in
        "Improved tangent estimate in the nudged elastic band method for finding
        minimum energy paths and saddle points"
        Graeme Henkelman and Hannes Jonsson
        J. Chem. Phys 113 (22), 9978 (2000)

        Parameters
        ----------
        
        central : float
            central image energy
        left : float
            left image energy
        right : float
            right image energy
        gleft : np.array
            gradient to left image (x_0 - x_left)
        gright : np.array
            gradient to right image (x_0 - x_right) 
        
        """

        vmax = max(abs(central - left), abs(central - right))
        vmin = min(abs(central - left), abs(central - right))

        tleft = gleft
        tright = -gright

        # special interpolation treatment for maxima/minima
        if (central >= left and central >= right) or (central <= left and central <= right):
            if left > right:
                t = vmax * tleft + vmin * tright
            else:
                t = vmin * tleft + vmax * tright
        # left is higher, take this one
        elif left > right:
            t = tleft
        # otherwise take right
        else:
            t = tright

        return t / norm(t)

    def NEBForce(self, isclimbing, image, left, right, greal, icenter):
        """
        Calculate NEB force for 1 image. That contains projected real force and spring force.

        The current implementation is the DNEB (doubly nudged elastic band) as described in

        "A doubly nudged elastic band method for finding transition states"
        Semen A. Trygubenko and David J. Wales
        J. Chem. Phys. 120, 2082 (2004); doi: 10.1063/1.1636455

        """

        # construct tangent vector
        d_left, g_left = self.distance(image[1], left[1], distance=True, grad=True)
        d_right, g_right = self.distance(image[1], right[1], distance=True, grad=True)
        self.distances[icenter - 1] = np.sqrt(d_left)
        self.distances[icenter] = np.sqrt(d_right)

        t = self.tangent(image[0], left[0], right[0], g_left, g_right)
        if isclimbing:
            return greal - 2. * np.dot(greal, t) * t

        if True:
            from . import _NEB_utils

            E, g_tot = _NEB_utils.neb_force(t, greal, d_left, g_left, d_right, g_right, self.k, self.dneb)
            if self.with_springenergy:
                return E, g_tot
            else:
                return 0., g_tot
        else:  # pragma: no cover

            # project out parallel part
            gperp = greal - np.dot(greal, t) * t

            # parallel part of spring force
            gs_par = self.k * (d_left - d_right) * t

            g_tot = gperp + gs_par

            if self.dneb:
                g_spring = self.k * (g_left + g_right)

                # perpendicular part of spring
                gs_perp = g_spring - np.dot(g_spring, t) * t
                # double nudging
                g_tot += gs_perp - np.dot(gs_perp, gperp) * gperp / np.dot(gperp, gperp)

            if self.with_springenergy:
                E = 0.5 / self.k * (d_left ** 2 + d_right ** 2)
            else:
                E = 0.
            return E, g_tot

    def _step(self, coords=None, energy=None, rms=None):
        self.step += 1
        if self.adjustk_freq <= 0:
            return
        if self.step % 5 == 0:
            self._adjust_k(coords)

    def _adjust_k(self, coords):
        tmp = self.coords.copy()
        tmp[1:self.nimages - 1, :] = coords.reshape(self.active.shape)

        d = []
        for i in range(0, self.nimages - 1):
            d.append(self.distance(tmp[i, :], tmp[i + 1, :], distance=True, grad=False)[0])

        d = np.array(np.sqrt(d))
        average_d = np.average(d)
        deviation = np.abs((d - average_d) / average_d)
        avdev = np.average(deviation)
        if avdev > self.adjustk_tol:
            self.k *= self.adjustk_factor
            if self.verbose >= 1:
                logger.debug("increasing DNEB force constant to %s", self.k)
        else:
            self.k /= self.adjustk_factor
            if self.verbose >= 1:
                logger.debug("decreasing DNEB force constant to %s", self.k)

    def MakeHighestImageClimbing(self):
        """
        Make the image with the highest energy a climbing image
        """
        emax = max(self.energies)
        for i in range(1, len(self.energies) - 1):
            if abs(self.energies[i] - emax) < 1e-10:
                self.isclimbing[i] = True

    def MakeAllMaximaClimbing(self):
        """
        Make all maxima along the neb climbing images.
        """
        for i in range(1, len(self.energies) - 1):
            if self.energies[i] > self.energies[i - 1] and self.energies[i] > self.energies[i + 1]:
                self.isclimbing[i] = True

    def printState(self):  # pragma: no cover
        """
        print the current state of the NEB.  Useful for bug testing
        """
        if self.printStateFile is None:
            # get a unique file name
            for n in range(500):
                self.printStateFile = "neb.EofS.%04d" % n
                if not os.path.isfile(self.printStateFile):
                    break
            logger.info("NEB energies will be written to %s", self.printStateFile)
        with open(self.printStateFile, "a") as fout:
            fout.write("\n\n")
            fout.write("#neb step: %d\n" % self.getEnergyCount)
            S = 0.
            for i in range(self.nimages):
                if i == 0:
                    dist = 0.
                else:
                    dist, t = self.distance(self.coords[i, :], self.coords[i - 1, :], grad=False)
                S += dist
                fout.write("%f %g\n" % (S, self.energies[i]))

    def copy(self):
        """ create a copy of the current neb """
        neb = copy.copy(self)
        neb.coords = neb.coords.copy()
        neb.energies = neb.energies.copy()
        neb.isclimbing = copy.deepcopy(neb.isclimbing)

        neb.active = neb.coords[1:neb.nimages - 1, :]
        return neb


#
# only testing stuff below here
#
#
# import nebtesting as test
#
# def nebtest(MyNEB=NEB, nimages=22):
# import pylab as pl
# from _interpolate import InterpolatedPath
# from pele.optimize import lbfgs_py
#    NEBquenchParams = dict()
#    NEBquenchParams["iprint"]=20
#    NEBquenchParams["debug"]=True
#    NEBquenchParams["maxErise"]=0.1
#    
#    x = np.arange(.5, 5., .05)
#    y = np.arange(.5, 5., .05)
#    z = np.zeros([len(x), len(y)])
#    potential = test.leps()
#    for i in range(0, len(x)):
#        for j in range(0, len(y)):
#                z[j, i] = potential.getEnergy([x[i], y[j]])
#    print "done"
#    #z[z>0.] = 0.
#    #pl.imshow(z)
#    #pl.show()
#    initial = np.array([.75, 2.]) #np.random.random(3)
#    final = np.array([2., .75]) #np.random.random(3)
##    from pele.optimize import quench
##    print "quench initial"
##    ret = quench.lbfgs_py(initial, potential.getEnergyGradient)
##    initial = ret[0]
##    print "quench final"
##    ret = quench.quench(final, potential.getEnergyGradient)
##    final = ret[0]
##    print "done with quenching"
##    print initial, final
#    #print "Initial: ", initial
#    #print "Final: ", final
#    #pl.imshow(z)
#
#    neb = MyNEB(InterpolatedPath(initial, final, nimages) ,potential, k=1000, dneb=False, quenchParams=NEBquenchParams)
#    tmp = neb.coords
#    energies_interpolate = neb.energies.copy()
#    pl.figure()
#    pl.subplot(2,2,1)
#    pl.subplots_adjust(wspace=0.3, left=0.05, right=0.95, bottom = 0.14)
#
#    pl.title("path")
#    #pl.contourf(x, y, z)
#    pl.pcolor(x, y, z, vmax=-0.5, cmap=pl.cm.PuBu)
#    pl.colorbar()
#    pl.plot(tmp[:, 0], tmp[:, 1], 'ko-')
#    print "optimizing NEB"
#    neb.optimize()#quenchRoutine=quench.fire)
#    print "done"
#    tmp = neb.coords
#    pl.plot(tmp[:, 0], tmp[:, 1], 'ro-')
#    pl.xlabel("x")
#    pl.ylabel("y")
#    pl.axis(xmin=0.5, xmax=2.5, ymin=0.5, ymax=2.5)
#
#    pl.subplot(1,2,2)
#    pl.title("energy")
#    pl.plot(energies_interpolate, 'ko-', label="interpolate")
#    pl.plot(neb.energies, 'ro-', label="neb")
#    pl.xlabel("image")
#    pl.ylabel("energy")
#    pl.legend(loc='best')
#    
#    pl.subplot(1,2,2)
#    pl.title("energy")
#    pl.plot(energies_interpolate, 'ko-', label="interpolate")
#    pl.plot(neb.energies, 'ro-', label="neb")
#    pl.xlabel("image")
#    pl.ylabel("energy")
#    pl.legend(loc='best')
#    pl.show()
#
#if __name__ == "__main__":
#    nebtest()

