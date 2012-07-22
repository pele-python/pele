import numpy as np
from pygmin.optimize import quench


def minima_from_ts(getEnergyGradient, xt, n=None, displace=1e-3, quenchRoutine=quench.mylbfgs):
    # if no direction is given, choose random direction
    if n==None:
        # TODO: replace by better algorithm with uniform sampling
        n = np.random.random(xt.shape)-0.5
    x1 = xt - displace*n
    x2 = xt + displace*n
    minimum1 = quenchRoutine(x1, getEnergyGradient)
    minimum2 = quenchRoutine(x2, getEnergyGradient)
    return minimum1, minimum2

    
    

if __name__ == '__main__':
    pass