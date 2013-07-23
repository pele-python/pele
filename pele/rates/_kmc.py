import random
import numpy as np

from pele.rates._rate_calculations import GraphReduction, graph_from_rates, _three_state_graph

def weighted_pick(weights):
    """sample uniformly from the objects with given weights
    
    Parameters
    ----------
    weights : dict
        a dictionary with objects as the keys and weights as the values
    
    Returns
    -------
    the object selected uniformly from the weights
    """
    if len(weights) == 0:
        raise ValueError("weights must not have zero length")
    r = np.random.uniform(0., sum(weights.itervalues()))
    s = 0.0
#    print r, len(weights)    
    for u, w in weights.iteritems():
        s += w
        if r < s: return u
    return u


class KineticMonteCarlo(object):
    """class to do kinetic Monte Carlo runs
    
    Parameters
    ----------
    graph : networkx.Graph()
        see GraphReduction for a description of how to attach transition
        probabilities and waiting times
    """
    def __init__(self, graph, debug=False):
        self.graph = graph
        self.debug = debug
    
    def next(self, u):
        udata = self.graph.node[u]
        
        rates = {u:udata["P"]}
        for x, v, uvdata in self.graph.edges(u, data=True):
            assert u == x
            kuv = uvdata[GraphReduction.Pkey(u, v)]
            rates[v] = kuv
        
        unext = weighted_pick(rates)
#        print "rates", u, ":", rates, "chosen", unext
        
        return unext, udata["tau"]
        
    
    def first_passage(self, a, B, niter=1000):
        """start at A and stop when you get to B, return the time elapsed
        
        Parameters
        ----------
        a : node
            the start node
        B : iterable
            the set of final states
        
        Returns
        -------
        total_time, niter
        """
        current_state = a
        B = set(B)
        if self.debug: path = [current_state]
        total_time = 0.
        niter = 0
        while current_state not in B:
            unext, time = self.next(current_state)
            current_state = unext
            if self.debug: path.append(current_state)
            total_time += time
            niter += 1
        
        if self.debug:
            print path[0],
            for u in path[1:]:
                print "->", u,
            print ""
        
        return total_time, niter
            

    def mean_first_passage_time(self, a, b, niter=1000):
        """compute the mean first passage time between nodes a and b
        """
        tavg = 0.
        for i in xrange(niter):
            time, count = self.first_passage(a, [b], niter)
    #        print time
            tavg += time
        
        tavg /= niter
        print "mean first passage time", tavg
        return tavg
    
#
# testing only below here
#

import unittest
class TestKMC(unittest.TestCase):
    def setUp(self):
        graph = _three_state_graph()
        self.kmc = KineticMonteCarlo(graph)
    
    def test_mfp(self):
        time = self.kmc.mean_first_passage_time(0,1,niter=1000)
        self.assertAlmostEqual(time, 1.0, 1)
        
        
def test():
    
        graph = _three_state_graph()
        kmc = KineticMonteCarlo(graph)
        kmc.mean_first_passage_time(0, 1, niter=10000)

if __name__ == "__main__":
    unittest.main()
