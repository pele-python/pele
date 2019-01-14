from __future__ import print_function
import numpy as np
import logging

from pele.optimize import mylbfgs
from pele.utils.rotations import vec_random_ndim

__all__ = ["minima_from_ts"]

logger = logging.getLogger("pele.connect")


def determine_pushoff(
        pot, coords, vec, stepmin=.01, **unused_kwargs):
    """apply the pushoff along the direction vec
    """
    if unused_kwargs:
        print("keywords:", list(unused_kwargs.keys()), "are obsolete and ignored in determine_pushoff")
    return coords + stepmin * vec / np.linalg.norm(vec)

#    if grad is None:
#        _, grad = pot.getEnergyGradient(coords)
#    gpar0 = np.dot(grad, vec) / vnorm
#    step = stepmin / vnorm
#    while True:
#        coords1 = step * vec + coords
#        _, grad = pot.getEnergyGradient(coords1)
#        rms = np.linalg.norm(grad) / np.sqrt(grad.size)
#        gpar = np.dot(grad, vec) / vnorm
#        finalstep = step
#        if verbose:
#            logger.debug("gpar %s %s %s %s", gpar, gpar0, abs((gpar - gpar0) / gpar0), step)
#        if abs((gpar - gpar0) / gpar0) > gdiff and rms > rms_min:
#            break
#        if step > stepmax:
#            logger.debug("warning: taking maximum step size from transition state %s", step)
#            break
#        step *= 2.
#    if verbose:
#        logger.debug("using pushoff of %s", finalstep)
#    return coords1


def minima_from_ts(pot, xt, n=None, quench=None, stepmin=0.01, **kwargs):
    """
    step off either side of a transition state and quench to find the minima
    
    Parameters
    ----------
    pot : potential object
    xt : array
        transition state coords
    n : array
        direction to step off
    quench : callable
        routine to use to do the quenching
    kwargs : dict
        parameters to pass to determine_pushoff
    """
    if n is None:
        # if no direction is given, choose random direction
        n = vec_random_ndim(xt.size)

    if quench is None:
        quench = lambda coords: mylbfgs(coords, pot)

    x1 = determine_pushoff(pot, xt, n, stepmin=stepmin, **kwargs)
    x2 = determine_pushoff(pot, xt, -n, stepmin=stepmin, **kwargs)
    minimum1 = quench(x1)
    minimum2 = quench(x2)
    return minimum1, minimum2


if __name__ == '__main__':
    pass

