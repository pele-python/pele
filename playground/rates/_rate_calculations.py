"""
routines to help with computing rates from one subset of a graph to another
"""

import networkx as nx
import numpy as np


class GraphReduction(object):
    """
    each transition state has data
    
    P : transition probability
    k : rate constant
    
    each node has data
    
    P : the probability to stay in this state
    tau : the occupation time
    """
    def __init__(self, graph):
        self.graph = graph
        
        self.debug = False
    
    def add_edge(self, u, v):
        """add an edge to the graph and initialize it with the appropriate data"""
        if self.debug: print "creating edge", u, v
        data = {self._Pkey(u, v):0., self._Pkey(v, u):0.}
        self.graph.add_edge(u, v, attr_dict=data)
        return self._get_edge_data(u, v)
    
    def _Pkey(self, u, v):
        return ("P", u, v)
    
    def _get_edge_data(self, u, v):
        try:
            return self.graph[u][v]
        except KeyError:
            return self.graph[v][u]

    
    def _update_edge(self, u, v, uxdata, vxdata, x, xdata):
        """
        update the probabilities of transition between u and v
        
        Puv -> Puv + Pux * Pxv / (1-Pxx)
        """
        assert u != v
        # if the edge doesn't exists, create it.
        try:
            uvdata = self._get_edge_data(u, v)
        except KeyError:
            uvdata = self.add_edge(u, v)
        
        
        Pux = uxdata[self._Pkey(u, x)]
        Pxu = uxdata[self._Pkey(x, u)]

        Pvx = vxdata[self._Pkey(v, x)]
        Pxv = vxdata[self._Pkey(x, v)]
        
        Pxx = xdata["P"]
        
        if self.debug:
            Puvold = uvdata[self._Pkey(u, v)]
            Pvuold = uvdata[self._Pkey(v, u)]
            
        uvdata[self._Pkey(u, v)] += Pux * Pxv / (1.-Pxx)
        uvdata[self._Pkey(v, u)] += Pvx * Pxu / (1.-Pxx)
        
        if self.debug:
            print "updating edge data", u, v, Puvold, "->", uvdata[self._Pkey(u, v)]
            print "updating edge data", v, u, Pvuold, "->", uvdata[self._Pkey(v, u)]
    
    def _update_node(self, u, x, xdata, uxdata):
        """
        update the waiting time and Puu for node u
        
        tauu -> tauu + Pux * taux / (1-Pxx)
        
        Puu -> Puu + Pux * Pxu / (1-Pxx)
        """
        taux = xdata["tau"]
        Pxx = xdata["P"]
        
#        uxdata = self._get_edge_data(u, x)
        udata = self.graph.node[u]

        Pxu = uxdata[self._Pkey(x, u)]
        Pux = uxdata[self._Pkey(u, x)]
        
        # update the waiting time at u
        udata["tau"] += Pux * taux / (1.-Pxx)
        
        # update Puu
        if self.debug:
            Pold = udata["P"]
        udata["P"] += Pux * Pxu / (1.-Pxx)
        
        if self.debug:
            print "updating node data", u, Pold, "->", udata["P"]
        
            
    
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
            print "  Pxv", uvdata[self._Pkey(u, v)], ": v =", v
#            assert 1 >= uvdata[self._Pkey(u, v)] >= 0
#            assert 1 >= uvdata[self._Pkey(v, u)] >= 0
            total_prob += uvdata[self._Pkey(u, v)]
        
        print "  total prob", total_prob
    
    def _check_node(self, u, verbose=True):
        udata = self.graph.node[u]  
#        print "checking node", u
        assert udata["tau"] > 0
        assert udata["P"] > 0
        
        total_prob = udata["P"]
        for x, v, uvdata in self.graph.edges(u, data=True):
            assert 1 >= uvdata[self._Pkey(u, v)] >= 0
            assert 1 >= uvdata[self._Pkey(v, u)] >= 0
            total_prob += uvdata[self._Pkey(u, v)]
        
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

def _make_random_graph(nnodes=16):
    L = int(np.sqrt(nnodes))
    graph = nx.grid_graph([L,L], periodic=False)
    print graph
    tool = GraphReduction(graph)
    
    for u, udata in graph.nodes(data=True):
        udata["P"] = np.random.rand()
        udata["tau"] = np.random.rand()
    
    for u, v, uvdata in graph.edges_iter(data=True):
        uvdata[tool._Pkey(u, v)] = np.random.rand()
        uvdata[tool._Pkey(v, u)] = np.random.rand()
    
    # now normalize
    for u, udata in graph.nodes(data=True):
        total_prob = udata["P"]
        edges = [(v, data) for x, v, data in graph.edges(u, data=True)]
        for v, uvdata in edges:
            assert v != u
            total_prob += uvdata[tool._Pkey(u, v)]
        
        udata["P"] /= total_prob
        for v, uvdata in edges:
            uvdata[tool._Pkey(u, v)] /= total_prob


    tool.check_graph()
    return graph

def test(nnodes=16):
    graph = _make_random_graph(nnodes)
    reducer = GraphReduction(graph)

    x = graph.nodes()[0]
    print "removing node", x
    reducer.remove_node(x)
    reducer.check_graph()
    
    x = graph.nodes()[0]
    print "removing node", x
    reducer.remove_node(x)
    reducer.check_graph()
    


if __name__ == "__main__":
    test()
    
    