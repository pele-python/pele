from __future__ import division
import numpy as np
from pele.potentials import BasePotential

class MLCost(BasePotential):
    """
    Cost function to be used for maximum likelihood optimization.
    
    Parameters
    ----------
    data : array of floats
        The observed data.
    log_probf : callable `log_probf(data, parameters)`
        Log of probability function which takes an array of data and an array
        of parameters and returns a float value
    probf : callable `probf(data, parameters)`
        Probability function which takes an array of data and an array
        of parameters and returns a float value
    
    Examples
    --------
    To get maximum likelihood estimates of the parameters, do e.g.:
    
        pot = MLCost(observed, probf=gauss)
        optimizer = LBFGS_CPP(parameters, pot)
        result = optimizer.run()
        opt_parameters = result.coords
    """
    def __init__(self, data, log_probf=None, probf=None):
        self.log_probf = log_probf
        self.probf = probf
        if self.log_probf == None and self.probf == None:
            raise Exception("provide either log_probf or probf")
        self.data = np.asarray(data)

    def getEnergy(self, parameters):
        if self.probf != None:
            return -np.sum(np.log(self.probf(self.data, np.asarray(parameters))))
        return -np.sum(self.log_probf(self.data, np.asarray(parameters)))
