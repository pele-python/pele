import numpy as np
from pygmin import defaults
import pygmin.optimize.quench as quench

__all__ = ["minima_from_ts"]

def determinePushoff(
            getEnergyGradient, coords, vec, gdiff=100., stepmin = 1e-5, stepmax = 0.2,
            grad = None):
    """
    determine a good stepsize to step off a transition state
    
    try to choose a stepsize so (gradpar-gradpar_ts) / gradpar_ts >= gdiff
        where gradpar and gradpar_ts is the parallel component of the gradient
        after the step and at the top of the transition state
    """
    vnorm = np.linalg.norm(vec)
    if grad is None:
        e, grad = getEnergyGradient(coords)
    gpar0 = np.dot(grad, vec)/vnorm
    step = stepmin / vnorm
    #print "gpar0", gpar0
    while True:
        coords1 = step * vec + coords
        e, grad = getEnergyGradient(coords1)
        gpar = np.dot(grad, vec)/vnorm
        #print "gpar", gpar, gpar0, abs((gpar-gpar0) / gpar0), step
        if abs((gpar-gpar0) / gpar0) > gdiff:
            break
        if step > stepmax:
            print "warning: taking maximum step size from transition state", step
            break
        step *= 2.
    return coords1

def minima_from_ts(getEnergyGradient, xt, n=None, displace=1e-3,
                   quenchRoutine=None, quenchParameters=dict()):
    # if no direction is given, choose random direction
    if n==None:
        # TODO: replace by better algorithm with uniform sampling
        n = np.random.random(xt.shape)-0.5
    
    quenchParameters = dict(defaults.quenchParams.items() + 
                            quenchParameters.items())
    if quenchRoutine==None:
        quenchRoutine = quench.mylbfgs
        #quenchRoutine = defaults.quenchRoutine
        #js850> this should be done more carefully
        quenchParameters = dict(quenchParameters.items() +
                                [("M",1)])
        
    
    #x1 = xt - displace*n
    x1 = determinePushoff(getEnergyGradient, xt, n)
    x2 = determinePushoff(getEnergyGradient, xt, -n)
    #x2 = xt + displace*n
    #e1,g1 = getEnergyGradient(x1)
    #e2,g2 = getEnergyGradient(x2)
    #print np.dot(g1,g2)
    minimum1 = quenchRoutine(x1, getEnergyGradient, **quenchParameters)
    minimum2 = quenchRoutine(x2, getEnergyGradient, **quenchParameters)
    return minimum1, minimum2

    
    

if __name__ == '__main__':
    pass
