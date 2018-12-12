from __future__ import division
import numpy as np
from pele.potentials import BasePotential
import copy


class MLCost(BasePotential):
    """
    Cost function to be used for maximum likelihood optimization.
    
    Parameters
    ----------self.assertLessEqual(confidence_intervalsl[i][0], par)
            self.assertLessEqual(par, confidence_intervalsl[i][1])
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
        if self.log_probf is None and self.probf is None:
            raise Exception("provide either log_probf or probf")
        self.data = np.asarray(data)

    def getEnergy(self, parameters):
        if self.probf is not None:
            return -np.sum(np.log(self.probf(self.data, np.asarray(parameters))))
        return -np.sum(self.log_probf(self.data, np.asarray(parameters)))

    def get_error_estimate(self, opt_parameters, log_l_variation=0.5):
        """
        Get confidence interval on ML estimated parameters.
        
        Parameters
        ----------
        opt_parameters : array of floats
            The optimum parameters as obtained from the ML method.
        log_l_variation : float, optional
            The variation in the log-likelihood upon changing the values of the
            parameters within the confidence interval:
            log_likelihood(interval[0]) <= log_likelihood(opt_parameters) - log_l_variation
            log_likelihood(interval[1]) <= log_likelihood(opt_parameters) - log_l_variation
        
        Examples
        --------
        To get an estimate on the confidence interval of the optimum parameters,
        do e.g.:
    
            error = pot.get_error_estimate(result)
        """
        self.opt_parameters = opt_parameters
        self.log_l_variation = log_l_variation
        self.minimum_cost = self.getEnergy(self.opt_parameters)
        self.cost_interval_edge = self.minimum_cost + self.log_l_variation
        return [self.get_interval(par_idx) for par_idx in range(len(self.opt_parameters))]

    def get_interval(self, par_idx):
        opt_par = self.opt_parameters[par_idx]
        step_size = opt_par / 1e5
        left_interval_edge = opt_par
        right_interval_edge = opt_par
        while True:
            trial_parameters = copy.copy(self.opt_parameters)
            trial_parameters[par_idx] = left_interval_edge
            if self.getEnergy(trial_parameters) < self.cost_interval_edge:
                left_interval_edge -= step_size
            else:
                break
        while True:
            trial_parameters = copy.copy(self.opt_parameters)
            trial_parameters[par_idx] = right_interval_edge
            if self.getEnergy(trial_parameters) < self.cost_interval_edge:
                right_interval_edge += step_size
            else:
                break
        return [left_interval_edge, right_interval_edge]
    
    
