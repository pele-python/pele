from __future__ import division
import numpy as np
from pele.potentials import BasePotential

class MLCost(BasePotential):
    """
    Cost function to be used for ML optimization.
    Takes probability function,
    its parameters,
    and observed data.
    To get maximum likelihood estimates of the parameters, do e.g.:
    pot = MLCost(gauss, observed)
    optimizer = LBFGS_CPP(parameters, pot)
    result = optimizer.run()
    opt_parameters = result.coords
    """
    def __init__(self, probf, data):
        self.probf = probf
        self.data = np.array(data)
    def getEnergy(self, coords):
        coords = np.array(coords)
        return -np.sum(np.log(self.probf(self.data, coords)))
