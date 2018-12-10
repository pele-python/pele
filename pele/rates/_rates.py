"""routines to compute rates from a database of minima
"""
from __future__ import print_function
import networkx as nx
import numpy as np

from pele.utils.disconnectivity_graph import database2graph
from pele.rates._ngt_cpp import NGT
from pele.rates._rates_linalg import reduce_rates, TwoStateRates, LinalgError

__all__ = ["RateCalculation", "RatesLinalg", "compute_committors"]

#def log_sum(log_terms):
#    """
#    return the log of the sum of terms whose logarithms are provided.
#    """
#    lmax = np.max(log_terms)
##    lsub = log_terms - lmax
#    result = np.log(np.sum(np.exp(log_terms - lmax))) + lmax
#    return result

def log_sum2(a, b):
    """
    return log( exp(a) + exp(b) )
    """
    if a > b:
        return a + np.log(1.0 + np.exp(-a + b) )
    else:
        return b + np.log(1.0 + np.exp(a - b) )


def log_equilibrium_occupation_probability(minimum, T):
    """return the log equilibrium occupation probability
    
    This is computed from the harmonic superposition approximation.  Some
    constants that are the same for all minima have been left out.
    """
    # warning, this has not been checked, there might be a bug
    return (-minimum.energy / T - np.log(minimum.pgorder)
            - 0.5 * minimum.fvib)


class _Minima2Rates(object):
    """prepare a list of transition states for a rate calculation
    
    Compute rate constants, and equilibrium occupation probabilities.
    Also, make sure the network of minima is sensible and all
    connected.
    """
    def __init__(self, transition_states, A, B, T=1., 
                  use_fvib=True):
        self.transition_states = transition_states
        self.A = set(A)
        self.B = set(B)
        self.beta = 1. / T
        self.use_fvib = use_fvib
        self.initialized = False

    def run(self):
        """compute rate constants from transition states and check and reduce graph
        
        also compute equilibrium occupation probabilities
        """
        self._initialized = True
        self._compute_rate_constants()
        self.rate_constants = reduce_rates(self.rate_constants, self.A, self.B)
        self._get_equilibrium_occupation_probabilities()


    def _get_local_log_rate(self, min1, min2, ts):
        """rate for going from min1 to min2
        
        TODO: should properly sum over all transition states between min1 and min2.
        
        from book Energy Landscapes page 387:
        
            sigma = 2 * min1.pgorder / ts.pgorder
            k(T) = sigma * exp((min1.fvib - ts.fvib)/2) * exp( -(ts.energy - min1.energy))
        
        where fvib is the log product of the harmonic mode frequencies  
        """
        if self.use_fvib:
            sigma = float(min1.pgorder) / (2. * np.pi * ts.pgorder)
            return (np.log(sigma) + (min1.fvib - ts.fvib)/2. 
                  - (ts.energy - min1.energy) * self.beta)
        else:
            return -(ts.energy - min1.energy) * self.beta
    
    def _transition_state_ok(self, ts):
        if ts.invalid:
            return False
        if ts.minimum1.invalid:
            return False
        if ts.minimum2.invalid:
            return False
        return True
    
    def _compute_rate_constants(self):
        """build the graph that will be used in the rate calculation"""
        # get rate constants over transition states.
        # sum contributions from multiple transition states between two minima.
        log_rates = dict()
        nts_skip_same = 0
        for ts in self.transition_states:
            if not self._transition_state_ok(ts):
                print("excluding invalid transition state from rate graph", ts.energy, ts.id()) 
                continue
            min1, min2 = ts.minimum1, ts.minimum2
            if min1.id() == 1664:
                print(min1.id(), min2.id())
            if min1 == min2:
#                print "skipping transition state with energy", ts.energy, "that connects minimum", min1.id(), "with itself", min2.id()
                nts_skip_same += 1
                continue
            u, v = min1, min2
            log_kuv = self._get_local_log_rate(min1, min2, ts)
            log_kvu = self._get_local_log_rate(min2, min1, ts)
            if (u,v) in log_rates:
#                print "found another transition state", u.id(), v.id()
                log_rates[(u,v)] = log_sum2(log_rates[(u,v)], log_kuv)
                log_rates[(v,u)] = log_sum2(log_rates[(v,u)], log_kvu)
            else:
                log_rates[(u,v)] = log_kuv
                log_rates[(v,u)] = log_kvu
        if nts_skip_same > 0:
            print("warning: not using", nts_skip_same, "transition states because they connect a minimum with itself")
            
        # should we remove the largest component from the rates to avoid underflow and overflow?
        # if so we need to multiply all the rates by this value
        if True:
            self.max_log_rate = max(log_rates.values())
            self.rate_norm = np.exp(-self.max_log_rate)
            print("time scale need to be multiplied by", 1./np.exp(self.max_log_rate))
            rates = dict(( (uv,np.exp(log_k - self.max_log_rate)) for uv, log_k in log_rates.items() ))
        else:
            self.max_log_rate = 0.
            self.rate_norm = 1.
            rates = dict(( (uv,np.exp(log_k)) for uv, log_k in log_rates.items() ))
        self.rate_constants = rates
        return self.rate_constants
        
    def _log_equilibrium_occupation_probability(self, minimum):
        """return the log equilibrium occupation probability
        
        This is computed from the harmonic superposition approximation.  Some
        constants that are the same for all minima have been left out.
        """
        # warning, this has not been checked, there might be a bug
        return (-self.beta * minimum.energy - np.log(minimum.pgorder)
                - 0.5 * minimum.fvib)

    def _get_equilibrium_occupation_probabilities(self, all=True):
#        if len(self.A) == 1 and len(self.B) == 1:
#            self.weights = None
#            return self.weights
        if all:
            nodes = set()
            for ts in self.transition_states:
                nodes.add(ts.minimum1)
                nodes.add(ts.minimum2)
        else:
            nodes = self.A.union(self.B)
        log_weights = dict()
        for m in nodes:
            if m.invalid:
                continue
            log_weights[m] = self._log_equilibrium_occupation_probability(m)
        
        # normalize the weights to avoid overflow or underflow when taking the exponential
        weight_max = max(log_weights.values())
        self.weights = dict(( (x, np.exp(w - weight_max)) 
                              for x, w in log_weights.items()
                            ))
        return self.weights

class RateCalculation(object):
    """compute transition rates from a database of minima and transition states
    
    Parameters 
    ----------
    transition_states: iteratable
        list (or iterator) of transition states that define the rate network
    A, B : iterables
        groups of minima specifying the reactant and product groups.  The rates
        returned will be the rate from A to B and vice versa.
    T : float
        temperature at which to do the the calculation.  Should be in units of
        energy, i.e. include the factor of k_B if necessary.
    """
    def __init__(self, transition_states, A, B, T=1., 
                  use_fvib=True):
        self.minima2rates = _Minima2Rates(transition_states, A, B, T=T,
                                         use_fvib=use_fvib)

    def compute_rates(self):
        """compute the rates from A to B and vice versa"""
        self.minima2rates.run()
        self.reducer = NGT(self.minima2rates.rate_constants, 
                           self.minima2rates.A, self.minima2rates.B, 
                           weights=self.minima2rates.weights)
        self.reducer.compute_rates()

    def compute_rates_and_committors(self):
        """compute the rates from A to B and vice versa"""
        self.minima2rates.run()
        self.reducer = NGT(self.minima2rates.rate_constants, 
                           self.minima2rates.A, self.minima2rates.B, 
                           weights=self.minima2rates.weights)
        self.reducer.compute_rates_and_committors()

    def get_rate_AB(self):
        return self.reducer.get_rate_AB() / self.minima2rates.rate_norm

    def get_rate_BA(self):
        return self.reducer.get_rate_BA() / self.minima2rates.rate_norm

    def get_rate_AB_SS(self):
        return self.reducer.get_rate_AB_SS() / self.minima2rates.rate_norm

    def get_rate_BA_SS(self):
        return self.reducer.get_rate_BA_SS() / self.minima2rates.rate_norm
    
    def get_committors(self):
        committors = self.reducer.get_committors()
        return committors

class RatesLinalg(object):
    """this class duplicates the behavior in RateCalculation, but with the linalg solver"""
    _initialized = False
    _times_computed = False
    def __init__(self, transition_states, A, B, T=1., 
                  use_fvib=True):
        self.minima2rates = _Minima2Rates(transition_states, A, B, T=T,
                                         use_fvib=use_fvib)

    def initialize(self):
        self.minima2rates.run()
        self.two_state_rates = TwoStateRates(self.minima2rates.rate_constants, 
                                             self.minima2rates.A, self.minima2rates.B, 
                                             weights=self.minima2rates.weights)

    def compute_rates(self):
        if not self._initialized:
            self.initialize()
        if not self._times_computed:
            self.two_state_rates.compute_rates()
        self._times_computed = True
        return self.two_state_rates.get_rate_AB() / self.minima2rates.rate_norm
    
    def get_mfptimes(self):
        if not self._initialized:
            self.initialize()
        if not self._times_computed:
            self.two_state_rates.compute_rates()
        times = self.two_state_rates.mfpt_computer.mfpt_dict
        
        self.mfpt_dict = dict(( (m, t * self.minima2rates.rate_norm) 
                                for m,t in times.items() ))
        for m in self.minima2rates.B:
            self.mfpt_dict[m] = 0.
        return self.mfpt_dict
    
    def compute_committors(self):
        if not self._initialized:
            self.initialize()
        self.two_state_rates.compute_committors()
        self.committors = self.two_state_rates.committor_dict
        for m in self.minima2rates.A:
            self.committors[m] = 0.
        for m in self.minima2rates.B:
            self.committors[m] = 1.
        return self.committors
    
    def get_committors(self):
        return self.committors

def compute_committors(transition_states, A, B, T=1.):
    """compute and return the committor probabilites between A and B
    
    First try the linear algebra method.  If it fails, then use NGT
    """
    try:
        rcalc = RatesLinalg(transition_states, A, B, T=T)
        rcalc.compute_committors()
        return rcalc.get_committors()
    except LinalgError:
        print("sparse linear algebra method failed.  Using NGT instead")
        rcalc = RateCalculation(transition_states, A, B, T=T)
        rcalc.compute_rates_and_committors()
        return rcalc.get_committors()
        

