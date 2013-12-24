"""
routines to help with computing rates from one subset of a graph to another
"""
import itertools
from collections import defaultdict

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
    sumk = defaultdict(lambda: 0.)
    
    # compute the sum of the outgoing rates for each node
    for edge, rate in rates.iteritems():
        u, v = edge
        sumk[u] += rate
    
    
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
        
        self._final_Pxx = dict()
        self._final_tau = dict()
        
        self.debug = debug
        self.initial_check_graph()
        self.check_graph()
    
    def _phase_one_remove_intermediates(self):
        intermediates = set(self.graph.nodes())
        intermediates.difference_update(self.A)
        intermediates.difference_update(self.B)
        intermediates = list(intermediates)
        # The calculation is faster if we remove the nodes with the least edges first
        intermediates.sort(key=lambda x: self.graph.degree(x))
        
        for x in intermediates:
            self.remove_node(x)
    
    def get_final_rates(self):
        kAB = sum(( (1.-self._final_Pxx[x]) / self._final_tau[x]
                     for x in self.A))
        kBA = sum(( (1.-self._final_Pxx[x]) / self._final_tau[x]
                     for x in self.B))
        return kAB, kBA
    
    def renormalize(self):
        self._phase_one_remove_intermediates()
        
        self._phase_two()
        
        self.rateAB, self.rateBA = self.get_final_rates()
        return self.rateAB, self.rateBA
            
#         while len(self.A) > 1:
#             x = self.A.pop()
#             self.remove_node(x)
#             
#         while len(self.B) > 1:
#             x = self.B.pop()
#             self.remove_node(x)
#         
#         u = iter(self.A).next()
#         v = iter(self.B).next()
#         self.rateAB = self._get_rate(u, v)
#         self.rateBA = self._get_rate(v, u)
#         
#         if self.debug:
#             print "rate ", u, "->", v, self.rateAB
#             print "rate ", v, "->", u, self.rateBA
#         return self.rateAB, self.rateBA

    def _phase_two_group(self, full_graph, group):
        for a in group:
            self.graph = full_graph.copy()
            Acopy = set(group)
            Acopy.remove(a)
            while len(Acopy) > 0:
                x = Acopy.pop()
                self.remove_node(x)
            
            adata = self.graph.node[a]
            self._final_Pxx[a] = adata["P"]
            self._final_tau[a] = adata["tau"]
        

    def _phase_two(self):
        """
        in this second phase we deal with starting and ending sets that have more
        than 1 element.  This follows the text above equation 19 in Wales 2009.
        This is called after all the intermediates have been decimated.
        
        
        for each element a in set A, compute tauF_a and PF_aa by decimating all 
        nodes in A except a.  Then the final rate from A to B
        
        kAB = (1/p_eq_A) sum_a PF_aB / tauF_a * p_eq_a
        
        where in the above PF_aB = 1-PF_aa and p_eq_a is the equilibrium
        probability to be in state a and p_eq_A = sum_a p_eq_a
        
        The inverse rate is, symmetrically,
        
        kBA = (1/p_eq_B) sum_b PF_bA / tauF_b * p_eq_b
        
        """
        full_graph = self.graph.copy()
        self._phase_two_group(full_graph, self.A)
        self._phase_two_group(full_graph, self.B)
        
        # restore the full graph
        self.graph = full_graph
            
        
#         a = iter(self.A).next()
#         b = iter(self.B).next()
#         abdata = self._get_edge_data(a, b)
#         abdata2 = self._get_edge_data(a, b, graph=full_graph) 
#         print "before", abdata, abdata2
#         abdata2[self.Pkey(a, b)] = 0.222
#         print "after ", abdata, abdata2
 
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
    
    def _get_edge_data(self, u, v, graph=None):
        if graph is None: graph = self.graph
        try:
            return graph[u][v]
        except KeyError:
            return graph[v][u]

    
    def _update_edge(self, u, v, uxdata, vxdata, x, xdata):
        """
        update the probabilities of transition between u and v upon removing node x
        
        Puv -> Puv + Pux * Pxv / (1-Pxx)
        Pvu -> Pvu + Pvx * Pxu / (1-Pxx)
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
        
        # in the paper, to avoid numerical errors DJW computes 
        # 1-Pxx as sum_j Pxj if Pxx > .99         
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

    def _check_A_B_connected(self, connected_components):
        ccset = [set(c) for c in connected_components]
        
        ca_intersections = [c.intersection(self.A) for c in ccset]
        cb_intersections = [c.intersection(self.B) for c in ccset]

        sizes = [len(ca) for ca in ca_intersections if len(ca) > 0]
        if len(sizes) != 1:
            assert len(sizes) != 0
            print "warning, the reactant set is not fully connected"
            print "   ", [c for c in ca_intersections if len(c) > 0]
            raise Exception("the reactant set is not fully connected")

        sizes = [len(cb) for cb in cb_intersections if len(cb) > 0]
        if len(sizes) != 1:
            assert len(sizes) != 0
            print "warning, the product set is not fully connected"
            print "   ", [c for c in cb_intersections if len(c) > 0]
            raise Exception("the product set is not fully connected")
        
        AB_connected = False
        for ca, cb in itertools.izip(ca_intersections, cb_intersections):
            if len(ca) > 0 and len(cb) > 0:
                AB_connected = True
                break
        if not AB_connected:
            raise Exception("product and reactant sets are not connected")
        
            
        return False

    def initial_check_graph(self):
        for a in self.A:
            if not self.graph.has_node(a):
                raise Exception("an element in the reactant set is not in the graph")
        for b in self.B:
            if not self.graph.has_node(b):
                raise Exception("an element in the product set is not in the graph")
        
        assert len(self.A.intersection(self.B)) == 0

        # check A and B are connected
        cc = nx.connected_components(self.graph)
        if len(cc) != 1:
            print "warning, graph is not fully connected.  There are", len(cc), "components"
            self._check_A_B_connected(cc)
            
#          for a, b in itertools.product(self.A, self.B):
#             try: a, b
            

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
#     unittest.main()
    test()
    
    