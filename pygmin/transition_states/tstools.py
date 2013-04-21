import numpy as np
import logging

from pygmin.optimize import mylbfgs

__all__ = ["minima_from_ts"]

logger = logging.getLogger("pygmin.connect")

def determinePushoff(
            pot, coords, vec, gdiff=100., stepmin = 1e-5, stepmax = 0.2, verbose=False,
            grad = None):
    """
    determine a good stepsize to step off a transition state
    
    try to choose a stepsize so (gradpar-gradpar_ts) / gradpar_ts >= gdiff
        where gradpar and gradpar_ts is the parallel component of the gradient
        after the step and at the top of the transition state
    """
    vnorm = np.linalg.norm(vec)
    if grad is None:
        e, grad = pot.getEnergyGradient(coords)
    gpar0 = np.dot(grad, vec)/vnorm
    step = stepmin / vnorm
    #print "gpar0", gpar0
    while True:
        coords1 = step * vec + coords
        e, grad = pot.getEnergyGradient(coords1)
        gpar = np.dot(grad, vec) / vnorm
        finalstep = step
        if verbose:
            logger.debug("gpar %s %s %s %s", gpar, gpar0, abs((gpar-gpar0) / gpar0), step)
        if abs((gpar-gpar0) / gpar0) > gdiff:
            break
        if step > stepmax:
            logger.debug("warning: taking maximum step size from transition state %s", step)
            break
        step *= 2.
    if verbose:
        logger.debug( "using pushoff of %s", finalstep)
    return coords1

def minima_from_ts(pot, xt, n=None,
                   quenchRoutine=None, quenchParams=dict(), **kwargs):
    """
    step off either side of a transition state and quench to find the minima
    
    Parameters
    ----------
    pot : potential object
    xt : array
        transition state coords
    stepmin : float
        initial guess for size of pushoff when stepping off the transition state
    stepmax : float
        maximum size of pushoff when stepping off the transition state
    gdiff : float
        criterion for choosing a step size.  Try to choose a stepsize so that:: 
        
            (gradpar - gradpar_ts) / gradpar_ts >= gdiff
        
        where gradpar and gradpar_ts is the parallel component of the gradient
        after the step and at the top of the transition state
    quenchRoutine : callable
        routine to use to do the quenching
    quenchParams : dict
        parameters to pass to quenchRoutine
    verbose : bool
    """
    # if no direction is given, choose random direction
    if n==None:
        # TODO: replace by better algorithm with uniform sampling
        n = np.random.random(xt.shape)-0.5
    
    quenchParams = dict(quenchParams.items())
    if quenchRoutine is None:
        quenchRoutine = mylbfgs
        #js850> this should be done more carefully
        quenchParams = dict(quenchParams.items() +
                                [("M",1)])
        
    
    #x1 = xt - displace*n
    x1 = determinePushoff(pot, xt, n, **kwargs)
    x2 = determinePushoff(pot, xt, -n, **kwargs)
    #x2 = xt + displace*n
    #e1,g1 = getEnergyGradient(x1)
    #e2,g2 = getEnergyGradient(x2)
    #print np.dot(g1,g2)
    minimum1 = quenchRoutine(x1, pot.getEnergyGradient, **quenchParams)
    minimum2 = quenchRoutine(x2, pot.getEnergyGradient, **quenchParams)
    return minimum1, minimum2

    
    

if __name__ == '__main__':
    pass
