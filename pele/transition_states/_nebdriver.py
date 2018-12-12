from __future__ import print_function
from __future__ import absolute_import
import logging
import numpy as np

from pele.transition_states import NEB
from pele.transition_states._NEB import distance_cart
from ._interpolate import InterpolatedPath, interpolate_linear
from pele.utils.events import Signal

__all__ = ["NEBDriver"]

logger = logging.getLogger("pele.connect.neb")


class NEBDriver(object):
    """ driver class for NEB

    The NEBDriver wraps calls for NEB from LocalConnect. The driver class is responsible for setting
    up the initial interpolation and optimizing the band.

    Parameters
    -----------
    potential :
        the potential object
    coords1, coords2 : array
        the structures to connect with the band
    k : float, optional
        the elastic band spring constant
    max_images : int
        the maximum number of NEB images
    image_density : float
        how many NEB images per unit distance to use.
    iter_density : float
        how many optimization iterations per unit distance to use.
    adjustk_freq : integer
        frequency to adjust k, set to zero to disable
    adjustk_tol : float
        tolerance for adjusting k up or down
    adjustk_factor : float
        the multiplicative factor used to adjust k
    dneb : bool
        use DNEB (Doubly-Nudged Elastic Band) rather than NEB
    reinterpolate : integer
        reinterpolate the path to achieve equidistant spacing every so many steps
    reinterpolate_tol : float
        tolerance for reinterpolation, only reinterpolate if relative change
        in nimages or distance variation are above tolerance
    adaptive_nimages : bool
        adjust number of images on reinterpolate to match image density
    adaptive_niter : bool
        adjust number of iterations if nimages is adjusted
    factor : float
        The number of images is multiplied by this factor.  If the number of
        images is already at it's maximum, then the number of iterations is
        multiplied by this factor instead
    verbose : integer
    interpolator : callable, optional
        the function used to do the path interpolation for the NEB
    NEBquenchParams : dict
        parameters passed to the minimizer
    kwargs : keyword options
        additional options are passed to the NEB class

    See Also
    ---------
    NEB
    InterpolatedPath

    """

    def __init__(self, potential, coords1, coords2,
                 k=100., max_images=50, image_density=10., iter_density=10.,
                 verbose=0, factor=1.05, NEBquenchParams=None, adjustk_freq=0,
                 adjustk_tol=0.1, adjustk_factor=1.05, dneb=True,
                 reinterpolate_tol=0.1,
                 reinterpolate=0, adaptive_nimages=False, adaptive_niter=False,
                 interpolator=interpolate_linear, distance=distance_cart, **kwargs):

        self.potential = potential
        self.interpolator = interpolator

        self.verbose = verbose
        self.distance = distance
        self.factor = factor
        self.max_images = max_images
        self.min_images = 3
        self.image_density = image_density
        self.iter_density = iter_density
        self.update_event = Signal()
        self.coords1 = coords1
        self.coords2 = coords2
        self.reinterpolate = reinterpolate
        self.reinterpolate_tol = reinterpolate_tol
        self._nebclass = NEB
        self._kwargs = kwargs.copy()
        self.k = k
        self.adaptive_images = adaptive_nimages
        self.adaptive_niter = adaptive_niter

        self._kwargs["adjustk_freq"] = adjustk_freq
        self._kwargs["adjustk_tol"] = adjustk_tol
        self._kwargs["adjustk_factor"] = adjustk_factor

        self._kwargs["dneb"] = dneb

        if NEBquenchParams is None:
            NEBquenchParams = {'tol': 1e-2, 'maxErise': 100.0}

        self.quenchParams = NEBquenchParams

        self.prepared = False

    @classmethod
    def params(cls, obj=None):
        if obj is None:
            obj = NEBDriver(None, None, None)
        ''' return the the parameters of current instance '''
        params = obj._kwargs.copy()
        params["k"] = obj.k
        params["image_density"] = obj.image_density
        params["max_images"] = obj.max_images
        params["iter_density"] = obj.iter_density
        params["reinterpolate"] = obj.reinterpolate
        params["reinterpolate_tol"] = obj.reinterpolate_tol
        params["adaptive_nimages"] = obj.adaptive_images
        params["adaptive_niter"] = obj.adaptive_niter

        params["verbose"] = obj.verbose
        params["NEBquenchParams"] = obj.quenchParams.copy()
        # factor is not here, todo, move this out of constuctor
        params["interpolator"] = obj.interpolator
        params["distance"] = obj.distance

        return params

    def prepare(self, path=None):
        self.prepared = True

        if path is None:
            path = self.generate_path(self.coords1, self.coords2)
        self.path = path

        self.nimages = len(self.path)
        self.steps_total = 0
        self.last_k = self.k

        energies = []
        for x in self.path:
            energies.append(self.potential.getEnergy(x))

        distances = []
        for i in range(len(self.path) - 1):
            distances.append(np.sqrt(self.distance(self.path[i], self.path[i + 1])[0]))

        self.update_event(path=np.array(self.path), energies=np.array(energies),
                          distances=np.array(distances), stepnum=self.steps_total,
                          rms=1.0, k=self.last_k, event="initial")

    def run(self):
        # determine the number of iterations                
        if not self.prepared:
            self.prepare()

        quenchParams = self.quenchParams.copy()
        # if nimages is already max_images then increasing the number
        # of images with factor will have no effect.  so double the number of steps instead
        niter = int(self.iter_density * self.nimages)
        if self.factor > 1. and self.nimages == self.max_images and self.max_images > 0:
            niter *= self.factor
            niter = int(niter)

        quenchParams["nsteps"] = niter

        if self.verbose >= 0:
            logger.info("    NEB: nimages   %s", self.nimages)
            logger.info("    NEB: nsteps    %s", niter)
            logger.info("    NEB: verbosity %s", self.verbose)

        if self.reinterpolate > 0:
            quenchParams["nsteps"] = min(self.reinterpolate, niter)

        self.niter = niter
        while True:
            # set up the NEB 
            k = self.last_k
            neb = self._nebclass(self.path, self.potential, k=k,
                                 quenchParams=quenchParams, verbose=self.verbose,
                                 distance=self.distance, **self._kwargs)
            self.neb = neb
            neb.events.append(self._process_event)

            # optimize the NEB            
            res = neb.optimize()
            self.last_k = neb.k

            # check if we're finished
            self.steps_total += res.nsteps
            if res.success or self.steps_total >= self.niter:
                res.nsteps = self.steps_total
                self.path = res.path

                if self.verbose >= 0:
                    logger.info("NEB finished after %d steps, rms %e" % (res.nsteps, res.rms))

                self._send_finish_event(res)
                return neb

            # get the distances between each of the images
            distances = []
            for i in range(len(res.path) - 1):
                distances.append(np.sqrt(self.distance(res.path[i], res.path[i + 1])[0]))

            # reinterplate the path based on the distances
            path = self._reinterpolate(res.path, distances)
            self.path = path

            # update the number of iterations if required
            if self.adaptive_niter:
                self.niter = int(self.iter_density * len(path))
                if self.factor > 1. and len(path) == self.max_images and self.max_images > 0:
                    self.niter = int(self.niter * self.factor)

            if self.verbose >= 1:
                logger.info("NEB reinterpolating path, %d images, niter is %d" % (len(path), self.niter))


    def generate_path(self, coords1, coords2):
        # determine the number of images to use
        dist, tmp = self.distance(coords1, coords2)
        dist = np.sqrt(dist)
        nimages = int(max(1., dist) * self.image_density * self.factor)
        if self.max_images > 0:
            nimages = min(nimages, self.max_images)
        if nimages < self.min_images:
            nimages = int(self.min_images)
        path = InterpolatedPath(coords1, coords2, nimages, interpolator=self.interpolator)

        return [x for x in path]

    def _reinterpolate(self, path, distances):
        average_d = np.average(distances)
        deviation = np.abs((distances - average_d) / average_d)
        avdev = np.average(deviation)

        acc_dist = np.sum(distances)
        nimages = len(path)
        if self.adaptive_images:
            nimages = int(int(max(1., acc_dist) * self.image_density * self.factor))
        if self.max_images > 0:
            nimages = int(min(nimages, self.max_images))
        if nimages < self.min_images:
            nimages = int(self.min_images)

        # only reinterpolate if above tolerance
        if (avdev < self.reinterpolate_tol and
                    abs(float(nimages - len(path)) / float(nimages)) < self.reinterpolate_tol):
            print("no reinterpolation needed")
            return path

        newpath = []
        newpath.append(path[0].copy())

        icur = 0
        s_cur = 0.
        s_next = distances[icur]

        for i in range(1, nimages - 1):
            s = float(i) * acc_dist / (nimages - 1)
            while s > s_next:
                icur += 1
                s_cur = s_next
                s_next += distances[icur]

            t = (s - s_cur) / (s_next - s_cur)
            newpath.append(self.interpolator(path[icur], path[icur + 1], t))
        newpath.append(path[-1].copy())
        return newpath

    def _process_event(self, path=None, energies=None, distances=None, stepnum=None, rms=None):
        self.update_event(path=path, energies=energies,
                          distances=distances, stepnum=stepnum + self.steps_total,
                          rms=rms, k=self.neb.k, event="update")

    def _send_finish_event(self, res):
        distances = []
        for i in range(len(res.path) - 1):
            distances.append(np.sqrt(self.distance(res.path[i], res.path[i + 1])[0]))

        self.update_event(path=res.path, energies=res.energy,
                          distances=np.array(distances), stepnum=res.nsteps,
                          rms=res.rms, k=self.neb.k, event="final")
        

