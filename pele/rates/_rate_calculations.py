"""
routines for computing rates from one subset of a graph to another
"""
import itertools
from collections import defaultdict

import networkx as nx
import numpy as np

def kmcgraph_from_rates(rates):
    """create a graph for input to GraphReduction from a dictionary of rates
    
    Parameters
    ----------
    rates : dict
        a dictionary of rates.  the keys are tuples of nodes (u,v), the values
        are the rates.
        
            rate_uv = rate[(u,v)]
    """
    graph = nx.DiGraph()
    sumk = defaultdict(lambda: 0.)
    
    # compute the sum of the outgoing rates for each node
    for edge, rate in rates.iteritems():
        u, v = edge
        sumk[u] += rate
    
    
    # add nodes to the rate graph and assign waiting time and Puu
    for u, sumk_u in sumk.iteritems():
        tau = 1. / sumk_u
        Puu = 0.
        graph.add_node(u, tau=tau)
        graph.add_edge(u, u, P=Puu)
    
    # add edges to rate graph and assign transition probabilities
    for edge, rate in rates.iteritems():
        u, v = edge
        tau_u = graph.node[u]["tau"]
        Puv =  rate * tau_u
        graph.add_edge(u, v, P=Puv)
    
    return graph
        
    
    
class GraphReduction(object):
    """
    class to apply the graph reduction method for finding transition rates between two groups of nodes
    
    Parameters
    ----------
    graph : networkx.DiGraph object
        A directed graph specifying the connectivity, the initial transition
        probabilities and the occupation times.  The graph must have all the
        data in the correct format.  Each node must have the following keys in
        their attributes dictionary::
            
            "tau" : occupation time
        
        Each edge between nodes u and v must have the following keys in their
        attributes dictionary:
        
            "P" : transition probability from u to v
        
    A, B : iterables
        Groups of nodes specifying the reactant and product groups.  The rates
        returned will be the rate from A to B and vice versa.
    weights : dict
        Dictionary with nodes as keys and weights as values.  The weights are
        the equilibrium occupation probabilities of the nodes in A and B.  They
        are used to do the weighted mean for the final average over inverse
        mean first passage times.
    
    Notes
    -----
    This follows the new graph transformation procedure (NGT) described by 
    David Wales, J. Chem. Phys., 2009 http://dx.doi.org/10.1063/1.3133782
    
    The rate, rAB computed by this calculation (returned by
    self.get_final_rates) is the inverse mean first passage time averaged over
    the states in A
    
    """
    def __init__(self, graph, A, B, debug=False, weights=None):
        self.graph = graph
        self.A = set(A)
        self.B = set(B)
        if weights is None:
            # everything has weight 1
            self.weights = defaultdict(lambda : 1.)
        else:
            self.weights = weights
        
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
        intermediates.sort(key=lambda x: self.graph.in_degree(x) + self.graph.out_degree(x))
        
        for x in intermediates:
            self.remove_node(x)
    
    def _get_final_rate(self, group):
        # should maybe be careful when Pxx is very close to 1.
        rate = sum(( (1.-self._final_Pxx[x]) / self._final_tau[x] * self.weights[x]
                     for x in group))
        norm = sum((self.weights[x] for x in group))
        return rate / norm
    
    def get_committor_probability(self, x):
        """return the committor probability for node x
        
        x must be in A or in B. If x is in A return the the probability
        that a trajectory starting at x gets to B before returning to x.
        If x is in B return the the probability
        that a trajectory starting at x gets to A before returning to x.
        """
        return 1. - self._final_Pxx[x]
    
    def get_rate_AB(self):
        """Return the transition rate from A to B
        
        This is the inverse mean first passage time averaged over the elements in A"""
        return self._get_final_rate(self.A)
    
    def get_rate_BA(self):
        """Return the transition rate from B to A
        
        This is the inverse mean first passage time averaged over the elements in B"""
        return self._get_final_rate(self.B)
    
    def compute_rates(self):
        """do the computation to compute the rates"""
        self._phase_one_remove_intermediates()
        
        self._phase_two()
        
#         self.rateAB, self.rateBA = self.get_final_rates()
#         return self.rateAB, self.rateBA

    def _phase_two_group(self, full_graph, group):
        """
        for each element a in the group, remove all other elements in
        the group then record the node attributes of a for later analysis.
        
        Note: This is a very inefficient way to do it.  If this becomes a 
        bottleneck it should be rewritten.
        """
        for a in group:
            self.graph = full_graph.copy()
            Acopy = set(group)
            Acopy.remove(a)
            while len(Acopy) > 0:
                x = Acopy.pop()
                self.remove_node(x)
            
            if self.graph.out_degree(a) <= 1:
                raise Exception("node %s is not connected" % (a))
            adata = self.graph.node[a]
            # in the paper, to avoid numerical errors DJW computes 
            # 1-Pxx as sum_j Pxj if Pxx > .99
            Paa = self._get_edge_data(a, a)["P"]
            if Paa > 0.999:
                print "warning, Pxx is very large (%s), numerical precision problems might be in your future" % (Paa)
            self._final_Pxx[a] = Paa

            self._final_tau[a] = adata["tau"]
        
    def _phase_two(self):
        """
        in this second phase we deal with reactant and product sets (A and B) that have more
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
            
    def _add_edge(self, u, v):
        """add an edge to the graph and initialize it with the appropriate data"""
        if self.debug: print "  creating edge", u, v
        self.graph.add_edge(u, v, P=0.)
        return self._get_edge_data(u, v)
    
    def _get_edge_data(self, u, v):
        return self.graph[u][v]

    def _update_edge(self, u, v, x, Pxx):
        """
        update the probabilities of transition between u and v upon removing node x
        
        Puv -> Puv + Pux * Pxv / (1-Pxx)
        
        u==v is fine, but u==x or v==x is not
        """
        assert u != x
        assert v != x
        try:
            uxdata = self._get_edge_data(u, x)
            xvdata = self._get_edge_data(x, v)
        except KeyError:
            # if either Pux or Pxv does not exist then nothing need be done
            return
        Pux = uxdata["P"]
        Pxv = xvdata["P"]
        # if the edge u, v doesn't exist, create it.
        try:
            uvdata = self._get_edge_data(u, v)
        except KeyError:
            uvdata = self._add_edge(u, v)
        
        if self.debug:
            Puvold = uvdata["P"]
        
        # update transition probability
        uvdata["P"] += Pux * Pxv / (1.-Pxx)
        
        if self.debug:
            print "  updating edge", u, "->", v, ":", Puvold, "->", uvdata["P"]
    
    def _update_node(self, u, x, tau_x, Pxx):
        """
        update the waiting time and Puu for node u upon removing node x
        
        tauu -> tauu + Pux * taux / (1-Pxx)
        
        Puu -> Puu + Pux * Pxu / (1-Pxx)
        """
        assert x != u
        
        udata = self.graph.node[u]

        Pux = self._get_edge_data(u, x)["P"]
        
        if self.debug:
            tauold = udata["tau"]

        # update the waiting time at u
        udata["tau"] += Pux * tau_x / (1.-Pxx)
        
        if self.debug:
            print "  updating node data", u, "tau", tauold, "->", udata["tau"]

    def remove_node(self, x):
        """
        remove node x from the graph and update the neighbors of x
        """
        neibs = set(self.graph.successors(x) + self.graph.predecessors(x))
        neibs.remove(x)
        tau_x = self.graph.node[x]["tau"]
        # in the paper, to avoid numerical errors DJW computes 
        # 1-Pxx as sum_j Pxj if Pxx > .99         
        Pxx = self._get_edge_data(x, x)["P"]
        if Pxx > 0.999:
            print "warning, Pxx is very large (%s), numerical precision problems might be in your future" % (Pxx)
        
        if self.debug:
            print "removing node", x, tau_x, Pxx

        # update node data
        for u in neibs:
            self._update_node(u, x, tau_x, Pxx)
        
        # update the edges between neighbors
        for u in neibs:
            for v in neibs:
                self._update_edge(u, v, x, Pxx)
        
        self.graph.remove_node(x)

    def _print_node_data(self, u):
        print "data from node x =", u
        udata = self.graph.node[u]  
#        print "checking node", u
        print "  taux",  udata["tau"]

        total_prob = 0.
        for x, v, uvdata in self.graph.out_edges_iter(u, data=True):
            Puv = uvdata["P"]
            print "  Pxv", Puv, ": v =", v
            total_prob += Puv
        
        print "  total prob", total_prob

    def _check_node(self, u, verbose=True):
        udata = self.graph.node[u]  
#        print "checking node", u
        assert udata["tau"] >= 0

        total_prob = 0.
        for x, v, uvdata in self.graph.out_edges_iter(u, data=True):
            Puv = uvdata["P"]
            assert 1 >= Puv >= 0
            total_prob += Puv

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

        # add node self loops with zero probability if they don't already exist
        for u in self.graph.nodes():
            try:
                self._get_edge_data(u, u)
            except:
                self._add_edge(u, u)
        
        assert len(self.A.intersection(self.B)) == 0

        # check A and B are connected
        undirected_graph = self.graph.to_undirected()
        cc = nx.connected_components(undirected_graph)
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

    
