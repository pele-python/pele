from __future__ import division
import numpy as np
from pele.potentials import BasePotential

class MLCost(BasePotential):
    """
    Cost function to be used for maximum likelihood optimization.
    
    Parameters
    ----------
    probf : callable `probf(data, parameters)`
        Probability function which takes an array of data and an array
        of parameters and returns a float value
    data : array of floats
        The observed data.
    
    Examples
    --------
    To get maximum likelihood estimates of the parameters, do e.g.:
    
        pot = MLCost(gauss, observed)
        optimizer = LBFGS_CPP(parameters, pot)
        result = optimizer.run()
        opt_parameters = result.coords
    """
    def __init__(self, probf, data):
        self.probf = probf
        self.data = np.array(data)

    def getEnergy(self, parameters):
        return -np.sum(np.log(self.probf(self.data, np.asarray(parameters))))
