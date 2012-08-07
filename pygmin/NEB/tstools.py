import numpy as np
from pygmin import defaults


def minima_from_ts(getEnergyGradient, xt, n=None, displace=1e-3,
                   quenchRoutine=None, quenchParameters=None):
    # if no direction is given, choose random direction
    if n==None:
        # TODO: replace by better algorithm with uniform sampling
        n = np.random.random(xt.shape)-0.5
    
    if quenchRoutine==None:
        quenchRoutine = defaults.quenchRoutine
    if quenchParameters==None:
        quenchParameters = defaults.quenchParams
        
    x1 = xt - displace*n
    x2 = xt + displace*n
    e1,g1 = getEnergyGradient(x1)
    e2,g2 = getEnergyGradient(x2)
    print np.dot(g1,g2)
    minimum1 = quenchRoutine(x1, getEnergyGradient, **quenchParameters)
    minimum2 = quenchRoutine(x2, getEnergyGradient, **quenchParameters)
    return minimum1, minimum2

    
    

if __name__ == '__main__':
    pass