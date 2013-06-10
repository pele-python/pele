import networkx as nx
import numpy as np

from pele.utils.disconnectivity_graph import database2graph
from _rate_calculations import GraphReduction

class RateCalculation(object):
    """compute transition rates from a database of minima and transition states
    
    Parameters 
    ----------
    database : pele Database object
    A, B : iterables
        groups of minima specifying the reactant and product groups.  The 
        rates returned will be the rate from A to B and vice versa.
    T : float
        temperature at which to do the the calculation.  Should be in units
        of energy, i.e. include the factor of k_B if necessary.
    """
    def __init__(self, database, A, B, T):
        self.database = database
        self.A = set(A)
        self.B = set(B)
        self.beta = 1. / T
        self.tsgraph = database2graph(database)
        
        
    def _reduce_tsgraph(self):
        """remove nodes from tsgraph that are not connected to A"""
        # get all nodes that are connected to A 
        u = iter(self.A).next()
        nodes = nx.single_source_shortest_path_length(self.tsgraph, u)
        nodes = set(nodes.iterkeys())
        nodes.add(u)
        
        # rebuild tsgraph with only those nodes
        self.tsgraph = self.tsgraph.subgraph(nodes)
        
        # check set A is in the graph
        Adiff = self.A.difference(nodes)
        if len(Adiff) > 0:
            print "warning: the reactant set A is not fully connected"
            self.A = self.A.difference_update(nodes)
        assert len(self.A) > 0

        # check set B is in the graph
        Bdiff = self.B.difference(nodes)
        if len(Bdiff) > 0:
            print "warning: the product set B is not fully connected or is not connected to A"
            self.B = self.B.difference_update(nodes)
        assert len(self.B) > 0
    
    def _get_local_rate(self, min1, min2, ts):
        """rate for going from min1 to min2"""
        return np.exp( -(ts.energy - min1.energy) * self.beta)
    
    def _min2node(self, minimum):
        return minimum._id
    
    def _make_rate_graph(self):
        """build the graph that will be used in the rate calculation"""
        self.rate_graph = nx.Graph()
        
        rates = dict()
        sumk = dict([(self._min2node(m), 0.) for m in self.tsgraph.nodes()])
        
        # get local transition rates
        for min1, min2, data in self.tsgraph.edges_iter(data=True):
            u = self._min2node(min1)
            v = self._min2node(min2)
            ts = data["ts"]
            kuv = self._get_local_rate(min1, min2, ts)
            kvu = self._get_local_rate(min2, min1, ts)
            rates[(u,v)] = kuv
            rates[(v,u)] = kvu
            sumk[u] += kuv
            sumk[v] += kvu
            
        # add nodes to rate graph and assign waiting time and Puu
        for u, sumku in sumk.iteritems():
            tau = 1. / sumku
            Puu = 0.
            data = {"P":Puu, "tau":tau}
            self.rate_graph.add_node(u, attr_dict=data)
        
        # add edges to rate graph and assign transition probabilities
        for min1, min2 in self.tsgraph.edges_iter():
            u = self._min2node(min1)
            v = self._min2node(min2)
            tauu = self.rate_graph.node[u]["tau"]
            tauv = self.rate_graph.node[v]["tau"]
            
            Puv = tauu * rates[(u,v)]
            Pvu = tauv * rates[(v,u)]
            
            data = {GraphReduction.Pkey(u, v):Puv,
                    GraphReduction.Pkey(v, u):Pvu}
            
            self.rate_graph.add_edge(u, v, attr_dict=data)

        self.Anodes = set([self._min2node(m) for m in self.A]) 
        self.Bnodes = set([self._min2node(m) for m in self.B]) 

    def compute_rates(self):
        """compute the rates from A to B and vice versa"""
        self._reduce_tsgraph()
        self._make_rate_graph()
        reducer = GraphReduction(self.rate_graph, self.Anodes, self.Bnodes)
        self.rateAB, self.rateBA = reducer.renormalize()
        return self.rateAB, self.rateBA

#
# only testing stuff below here
#

def test():
    from pele.systems import LJCluster
    from pele.landscape import ConnectManager
    system = LJCluster(13)
    db = system.create_database("lj13.db")
    
    if False:
        bh = system.get_basinhopping(db, outstream=None)
        bh.run(50)
        
        manager = ConnectManager(db, strategy="gmin")
        for i in range(10):
            min1, min2 = manager.get_connect_job("gmin")
            connect = system.get_double_ended_connect(min1, min2, db, verbosity=0)
            connect.connect()
    
    A = [db.minima()[0]]
    B = [db.minima()[-1]]
    
    rcalc = RateCalculation(db, A, B, T=1.)
    rcalc.compute_rates()

if __name__ == "__main__":
    test()