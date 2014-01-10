"""routines to compute rates from a database of minima
"""
import networkx as nx
import numpy as np

from pele.utils.disconnectivity_graph import database2graph
from pele.rates._rate_calculations import GraphReduction, kmcgraph_from_rates

__all__ = ["RateCalculation"]

def log_sum(log_terms):
    """
    return the log of the sum of terms whose logarithms are provided.
    """
    lmax = np.max(log_terms)
#    lsub = log_terms - lmax
    result = np.log(np.sum(np.exp(log_terms - lmax))) + lmax
    return result

def log_sum2(a, b):
    """
    return log( exp(a) + exp(b) )
    """
    if a > b:
        return a + np.log(1.0 + np.exp(-a + b) )
    else:
        return b + np.log(1.0 + np.exp(a - b) )



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
    ndof : int
        number of vibrational degrees freedom.  Must be included if A or B contains
        more than 1 minimum.
    """
    def __init__(self, transition_states, A, B, T=1., ndof=None, 
                  use_fvib=True):
        self.transition_states = transition_states
        self.A = set(A)
        self.B = set(B)
        self.beta = 1. / T
        self.ndof = ndof
        self.use_fvib = use_fvib

        if self.ndof is None:
            if len(self.A) > 1 or len(self.B) > 1:
                raise ValueError("if A or B has more than 1 minimum you must pass ndof")

#    def _reduce_graph(self, graph):
#        """remove nodes from tsgraph that are not connected to A"""
#        # remove nodes with invalid minima
#        nodes = [m for m in self.graph.nodes() if m.invalid]
#        if len(nodes) > 0:
#            print "removing %s invalid minima from rate graph" % (len(nodes))
#            graph.remove_nodes_from(nodes)
#        
#        # remove invalid transition states
#        edges = [(u,v) for u,v,data in graph.edges_iter(data=True) if data["ts"].invalid]
#        if len(edges) > 0:
#            print "removing %s invalid transition states from rate graph" % (len(nodes))
#            graph.remove_edges_from(edges)
#        
#        
#        # get all nodes that are connected to A 
#        u = iter(self.A).next()
#        nodes = nx.node_connected_component(graph, u)
#        nodes = set(nodes)
#        nodes.add(u)
#
#        # rebuild tsgraph with only those nodes
#        graph = graph.subgraph(nodes)
#        
#        # check set A is in the graph
#        Adiff = self.A.difference(nodes)
#        if len(Adiff) > 0:
#            print "warning: the reactant set A is not fully connected"
#            self.A = self.A.difference_update(nodes)
#        assert len(self.A) > 0
#
#        # check set B is in the graph
#        Bdiff = self.B.difference(nodes)
#        if len(Bdiff) > 0:
#            print "warning: the product set B is not fully connected or is not connected to A"
#            self.B = self.B.difference_update(nodes)
#        if len(self.B) == 0:
#            raise Exception("the product set B is empty or not connected to A")
        
    
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
    
    def _min2node(self, minimum):
        return minimum._id
    
    def _transition_state_ok(self, ts):
        if ts.invalid:
            return False
        if ts.minimum1.invalid:
            return False
        if ts.minimum2.invalid:
            return False
        return True
        
    def _make_kmc_graph(self):
        """build the graph that will be used in the rate calculation"""
        # get rate constants over transition states.
        # sum contributions from multiple transition states between two minima.
        log_rates = dict()
        for ts in self.transition_states:
            if not self._transition_state_ok(ts):
                print "excluding invalid transition state from rate graph", ts.energy, ts._id 
                continue
            min1, min2 = ts.minimum1, ts.minimum2
            u = self._min2node(min1)
            v = self._min2node(min2)
            log_kuv = self._get_local_log_rate(min1, min2, ts)
            log_kvu = self._get_local_log_rate(min2, min1, ts)
            if (u,v) in log_rates:
                log_rates[(u,v)] = log_sum2(log_rates[(u,v)], log_kuv)
                log_rates[(v,u)] = log_sum2(log_rates[(v,u)], log_kvu)
            else:
                log_rates[(u,v)] = log_kuv
                log_rates[(v,u)] = log_kvu
            
        # should we remove the largest component from the rates to avoid underflow and overflow?
        # if so we need to multiply all the rates by this value
#        max_log_rate = max(logrates.itervalues())
#        print "time scale need to be multiplied by", np.exp(max_log_rate)
#        rates = dict(( (uv,np.exp(log_k - max_log_rate)) for uv, log_k in logrates.iteritems() ))
        rates = dict(( (uv,np.exp(log_k)) for uv, log_k in log_rates.iteritems() ))
        
        # make the rate graph from the rate constants
        self.kmc_graph = kmcgraph_from_rates(rates)

        # translate the product and reactant set into the new node definition  
        self.Anodes = set([self._min2node(m) for m in self.A]) 
        self.Bnodes = set([self._min2node(m) for m in self.B]) 

    def _log_equilibrium_occupation_probability(self, minimum):
        """return the log equilibrium occupation probability
        
        This is computed from the harmonic superposition approximation.  Some
        constants that are the same for all minima have been left out.
        """
        # warning, this has not been checked, there might be a bug
        return (-self.beta * minimum.energy - np.log(minimum.pgorder)
                - 0.5 * minimum.fvib - self.ndof * np.log(self.beta))

    def _get_equilibrium_occupation_probabilities(self):
        if self.ndof is None:
            assert len(self.A) == 1
            assert len(self.B) == 1
            self.weights = None
            return self.weights
        log_weights = dict()
        for m in self.A.union(self.B):
            x = self._min2node(m)
            log_weights[x] = self._log_equilibrium_occupation_probability(m)
        
        # normalize the weights to avoid overflow or underflow when taking the exponential
        weight_max = max(log_weights.itervalues())
        self.weights = dict(( (x, np.exp(w - weight_max)) 
                              for x, w in log_weights.iteritems()
                            ))
        return self.weights

    def compute_rates(self):
        """compute the rates from A to B and vice versa"""
#        self._reduce_tsgraph()
        self._make_kmc_graph()
        weights = self._get_equilibrium_occupation_probabilities()
        self.reducer = GraphReduction(self.kmc_graph, self.Anodes, self.Bnodes, weights=weights)
        self.reducer.compute_rates()
        self.rateAB = self.reducer.get_rate_AB()
        self.rateBA = self.reducer.get_rate_BA()
        return self.rateAB, self.rateBA

#
# only testing stuff below here
#

def test():
    from pele.systems import LJCluster
    from pele.landscape import ConnectManager
    from pele.thermodynamics import get_thermodynamic_information
    system = LJCluster(13)
    system.params.structural_quench_params.tol = 1e-6
    db = system.create_database("lj13.db")
    
    if db.number_of_minima() < 10:
        bh = system.get_basinhopping(db, outstream=None)
        bh.run(50)
        
        manager = ConnectManager(db, strategy="gmin")
        for i in range(10):
            min1, min2 = manager.get_connect_job("gmin")
            connect = system.get_double_ended_connect(min1, min2, db, verbosity=0)
            connect.connect()
        

    connect = system.get_double_ended_connect(db.minima()[0], db.minima()[-1], db, verbosity=0)
    connect.connect()
        
    
    A = [db.minima()[0]]
    B = [db.minima()[-1]]
    
    get_thermodynamic_information(system, db, nproc=2)
    
    graph = database2graph(db)
    rcalc = RateCalculation(graph, A, B, T=1.)
    rAB, rBA = rcalc.compute_rates()
    print "rates", rAB, rBA

if __name__ == "__main__":
    test()
