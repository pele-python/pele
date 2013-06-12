"""
routines to help with computing rates from one subset of a graph to another
"""

import networkx as nx
import numpy as np

def graph_from_rates(rates):
    """create a graph for input to GraphReduction from a dictionary of rates
    
    Parameters
    ----------
    rates : dict
        a dictionary of rates.  the keys are tuples of nodes (u,v), the
        values are the rates.
        
            rate_uv = rate[(u,v)]
    """
    graph = nx.Graph()
    sumk = dict()
    
    # compute the sum of the outgoing rates for each node
    for edge, rate in rates.iteritems():
        u, v = edge
        try:
            sumk[u] += rate
        except KeyError:
            sumk[u] = rate
    
    
    # add nodes to the rate graph and assign waiting time and Puu
    for u, sumk_u in sumk.iteritems():
        tau = 1. / sumk_u
        Puu = 0.
        data = {"P":Puu, "tau":tau}
        graph.add_node(u, attr_dict=data)
    
    # add edges to rate graph and assign transition probabilities
    for edge, rate in rates.iteritems():
        u, v = edge
        tau_u = graph.node[u]["tau"]
        Puv =  rate * tau_u
        data = {GraphReduction.Pkey(u, v):Puv}
        graph.add_edge(u, v, attr_dict=data)
    
    return graph
        
    
    
        
        


class GraphReduction(object):
    """
    class to apply the graph reduction method for finding transition rates between two groups of nodes
    
    Parameters
    ----------
    graph : networkx.Graph object
        an undirected graph specifying the connectivity, the initial transition probabilities
        and the occupation times.  The graph must have all the data in the correct format.
        Each node must have the following keys in their attributes dictionary::
            
            "P"   : the probability to stay in this state
            "tau" : occupation time
        
        Each edge between nodes u and v must have the following keys in their 
        attributes dictionary:
        
            ("P", u, v) : transition probability from u to v
            ("P", v, u) : transition probability from v to u
        
    A, B : iterables
        groups of nodes specifying the reactant and product groups.  The rates returned will be 
        the rate from A to B and vice versa.
    """
    def __init__(self, graph, A, B, debug=False):
        self.graph = graph
        self.A = set(A)
        self.B = set(B)
        
        self.debug = debug
        self.check_graph()
    
    def renormalize(self):
        intermediates = set(self.graph.nodes())
        intermediates.difference_update(self.A)
        intermediates.difference_update(self.B)
        
        for x in intermediates:
            self.remove_node(x)
            
        while len(self.A) > 1:
            x = self.A.pop()
            self.remove_node(x)
            
        while len(self.B) > 1:
            x = self.B.pop()
            self.remove_node(x)
        
        u = iter(self.A).next()
        v = iter(self.B).next()
        self.rateAB = self._get_rate(u, v)
        self.rateBA = self._get_rate(v, u)
        
        if self.debug:
            print "rate ", u, "->", v, self.rateAB
            print "rate ", v, "->", u, self.rateBA
        return self.rateAB, self.rateBA
            
    def _get_rate(self, u, v):
        uvdata = self._get_edge_data(u, v)
        Puv = uvdata[self.Pkey(u, v)]
        tauu = self.graph.node[u]["tau"]
        
        return Puv / tauu
        
    
    def _add_edge(self, u, v):
        """add an edge to the graph and initialize it with the appropriate data"""
        if self.debug: print "creating edge", u, v
        data = {self.Pkey(u, v):0., self.Pkey(v, u):0.}
        self.graph.add_edge(u, v, attr_dict=data)
        return self._get_edge_data(u, v)
    
    @staticmethod
    def Pkey(u, v):
        """return the edge attribute key for the transition probability from u to v"""
        return ("P", u, v)
    
    def _get_edge_data(self, u, v):
        try:
            return self.graph[u][v]
        except KeyError:
            return self.graph[v][u]

    
    def _update_edge(self, u, v, uxdata, vxdata, x, xdata):
        """
        update the probabilities of transition between u and v upon removing node x
        
        Puv -> Puv + Pux * Pxv / (1-Pxx)
        """
        assert u != v
        # if the edge doesn't exists, create it.
        try:
            uvdata = self._get_edge_data(u, v)
        except KeyError:
            uvdata = self._add_edge(u, v)
        
        
        Pux = uxdata[self.Pkey(u, x)]
        Pxu = uxdata[self.Pkey(x, u)]

        Pvx = vxdata[self.Pkey(v, x)]
        Pxv = vxdata[self.Pkey(x, v)]
        
        Pxx = xdata["P"]
        
        if self.debug:
            Puvold = uvdata[self.Pkey(u, v)]
            Pvuold = uvdata[self.Pkey(v, u)]
            
        uvdata[self.Pkey(u, v)] += Pux * Pxv / (1.-Pxx)
        uvdata[self.Pkey(v, u)] += Pvx * Pxu / (1.-Pxx)
        
        if self.debug:
            print "updating edge data", u, v, Puvold, "->", uvdata[self.Pkey(u, v)]
            print "updating edge data", v, u, Pvuold, "->", uvdata[self.Pkey(v, u)]
    
    def _update_node(self, u, x, xdata, uxdata):
        """
        update the waiting time and Puu for node u upon removing node x
        
        tauu -> tauu + Pux * taux / (1-Pxx)
        
        Puu -> Puu + Pux * Pxu / (1-Pxx)
        """
        taux = xdata["tau"]
        Pxx = xdata["P"]
        
#        uxdata = self._get_edge_data(u, x)
        udata = self.graph.node[u]

        Pxu = uxdata[self.Pkey(x, u)]
        Pux = uxdata[self.Pkey(u, x)]
        
        if self.debug:
            tauold = udata["tau"]
            Pold = udata["P"]

        # update the waiting time at u
        udata["tau"] += Pux * taux / (1.-Pxx)
        
        # update Puu
        udata["P"] += Pux * Pxu / (1.-Pxx)
        
        if self.debug:
            print "updating node data", u, "P", Pold, "->", udata["P"], "tau", tauold, "->", udata["tau"]
        
            
    
    def remove_node(self, x):
        """
        remove node x from the graph and update the neighbors of x
        """
        neibs = [(v, data) for a, v, data in self.graph.edges(x, data=True)]
        xdata = self.graph.node[x]
        
        # update node data
        for u, uxdata in neibs:
            assert u != x
            self._update_node(u, x, xdata, uxdata)
        
        # update the edges between neighbors
        for i in xrange(len(neibs)):
            u , udata = neibs[i]
            for j in xrange(i+1, len(neibs)):
                v, vdata = neibs[j]
                self._update_edge(u, v, udata, vdata, x, xdata)
        
        self.graph.remove_node(x)
    
    def _print_node_data(self, u):
        print "data from node x =", u
        udata = self.graph.node[u]  
#        print "checking node", u
        print "  taux",  udata["tau"]
        print "  Pxx",  udata["P"]
        
        total_prob = udata["P"]
        for x, v, uvdata in self.graph.edges(u, data=True):
            print "  Pxv", uvdata[self.Pkey(u, v)], ": v =", v
#            assert 1 >= uvdata[self.Pkey(u, v)] >= 0
#            assert 1 >= uvdata[self.Pkey(v, u)] >= 0
            total_prob += uvdata[self.Pkey(u, v)]
        
        print "  total prob", total_prob
    
    def _check_node(self, u, verbose=True):
        udata = self.graph.node[u]  
#        print "checking node", u
        assert udata["tau"] >= 0
        assert 1 >= udata["P"] >= 0
        
        total_prob = udata["P"]
        for x, v, uvdata in self.graph.edges(u, data=True):
            assert 1 >= uvdata[self.Pkey(u, v)] >= 0
            assert 1 >= uvdata[self.Pkey(v, u)] >= 0
            total_prob += uvdata[self.Pkey(u, v)]
        
        assert np.abs(total_prob - 1.) < 1e-6, "%s: total_prob %g" % (str(u), total_prob)

    def check_graph(self):
        for u in self.graph.nodes():
            try:
                self._check_node(u, verbose=False)
            except:
                self._print_node_data(u)
                raise

#
# only testing stuff below here
#

import unittest
class TestGraphReduction(unittest.TestCase):
    def setUp(self):
        tmatrix = [ [0., 1., 1.,], [1., 0., 1.,], [1., 1., 0.] ]
        rates = dict()
        for i in range(3):
            for j in range(3):
                if i != j:
                    rates[(i,j)] = tmatrix[i][j]
    
        self.graph = graph_from_rates(rates)

    def test_rate(self):
        reducer = GraphReduction(self.graph, [0], [1], debug=True)
        reducer.check_graph()
        rAB, rBA = reducer.renormalize()
        print rAB, rBA
        reducer.check_graph()
        print reducer.graph.number_of_nodes(), reducer.graph.number_of_edges()

        

def _make_random_graph(nnodes=16):
    L = int(np.sqrt(nnodes))
    graph = nx.grid_graph([L,L], periodic=False)
    print graph
    tool = GraphReduction
    
    for u, udata in graph.nodes(data=True):
        udata["P"] = np.random.rand()
        udata["tau"] = np.random.rand()
    
    for u, v, uvdata in graph.edges_iter(data=True):
        uvdata[tool.Pkey(u, v)] = np.random.rand()
        uvdata[tool.Pkey(v, u)] = np.random.rand()
    
    # now normalize
    for u, udata in graph.nodes(data=True):
        total_prob = udata["P"]
        edges = [(v, data) for x, v, data in graph.edges(u, data=True)]
        for v, uvdata in edges:
            assert v != u
            total_prob += uvdata[tool.Pkey(u, v)]
        
        udata["P"] /= total_prob
        for v, uvdata in edges:
            uvdata[tool.Pkey(u, v)] /= total_prob


#    tool.check_graph()
    return graph




def test(nnodes=36):
    graph = _make_random_graph(nnodes)
    reducer = GraphReduction(graph, [], [])

    x = graph.nodes()[0]
    print "removing node", x
    reducer.remove_node(x)
    reducer.check_graph()
    
    x = graph.nodes()[0]
    print "removing node", x
    reducer.remove_node(x)
    reducer.check_graph()
    
    A = set(graph.nodes()[:2])
    B = set(graph.nodes()[-2:])
    print "A B", A, B
    reducer = GraphReduction(graph, A, B)
    rAB, rBA = reducer.renormalize()
    reducer.check_graph()
    print "number of nodes", graph.number_of_nodes()
    print "rates", rAB, rBA

if __name__ == "__main__":
    unittest.main()
#    test()
    
    