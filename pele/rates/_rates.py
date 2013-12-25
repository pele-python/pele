"""routines to compute rates from a database of minima
"""
import networkx as nx
import numpy as np

from pele.utils.disconnectivity_graph import database2graph
from pele.rates._rate_calculations import GraphReduction, graph_from_rates

__all__ = ["RateCalculation"]

class RateCalculation(object):
    """compute transition rates from a database of minima and transition states
    
    Parameters 
    ----------
    graph : networkx.Grapph
        transition state graph with the minima as nodes and the transition states
        attached to the edges. You can use the function database2graph() to 
        create this from a database.:
        
            from pele.utils.disconnectivity_graph import database2graph
            graph = database2graph(database)

    A, B : iterables
        groups of minima specifying the reactant and product groups.  The 
        rates returned will be the rate from A to B and vice versa.
    T : float
        temperature at which to do the the calculation.  Should be in units
        of energy, i.e. include the factor of k_B if necessary.
    """
    def __init__(self, graph, A, B, T, use_fvib=True):
        self.tsgraph = graph
        self.A = set(A)
        self.B = set(B)
        self.beta = 1. / T
        self.use_fvib = use_fvib
        
    def _reduce_tsgraph(self):
        """remove nodes from tsgraph that are not connected to A"""
        # get all nodes that are connected to A 
        u = iter(self.A).next()
        nodes = nx.node_connected_component(self.tsgraph, u)
        nodes = set(nodes)
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
        
#         if len(self.A) > 1 or len(self.B) > 1:
#             print "warning rates between minima in the product and reactant set should be set to zero but this isn't implemented yet"
#             # TODO : implement this
    
    def _get_local_rate(self, min1, min2, ts):
        """rate for going from min1 to min2
        
        TODO: should properly sum over all transition states between min1 and min2.
        
        from book Energy Landscapes page 387:
        
            sigma = 2 * min1.pgorder / ts.pgorder
            k(T) = sigma * exp(min1.fvib - ts.fvib) * exp( -(ts.energy - min1.energy))
        
        where fvib is the log product of the harmonic mode frequencies  
        """
        if self.use_fvib:
            sigma = 2. * min1.pgorder / ts.pgorder
            return sigma * np.exp(min1.fvib - ts.fvib 
                                  - (ts.energy - min1.energy) * self.beta )
        else:
            return np.exp( -(ts.energy - min1.energy) * self.beta)
    
    def _min2node(self, minimum):
        return minimum._id
    
    def _make_rate_graph(self):
        """build the graph that will be used in the rate calculation"""
        # get rate constants over transition states
        rates = dict()
        for min1, min2, data in self.tsgraph.edges_iter(data=True):
            u = self._min2node(min1)
            v = self._min2node(min2)
            ts = data["ts"]
            kuv = self._get_local_rate(min1, min2, ts)
            kvu = self._get_local_rate(min2, min1, ts)
            rates[(u,v)] = kuv
            rates[(v,u)] = kvu
        
        # make the rate graph from the rate constants
        self.rate_graph = graph_from_rates(rates)

        # translate the product and reactant set into the new node definition  
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
