'''Wrapper to represent a storage class as a graph'''
import networkx as nx

__all__ = ["TSGraph", "Graph"]


class _ConnectedComponents():
    """
    a class to build and maintain a dict of connected components
    
    This allows connections to be determined much more rapidly because
    the breadth first search algorithm (connected_components, bidirectional_dijkstra) 
    need be run only infrequently.
    """
    def __init__(self, graph):
        self.graph = graph
        self.need_rebuild = True
    
    def build(self):
        cc = nx.connected_components(self.graph)
        component = dict()
        #cc is a list of lists
        for i in range(len(cc)):
            nodes = cc[i]
            for n in nodes:
                component[n] = i
        self.component = component
        self.need_rebuild = False
    
    def setRebuild(self):
        self.need_rebuild = True

    def areConnected(self, min1, min2):
        if self.need_rebuild:
            self.build()
        return self.component[min1] == self.component[min2]
            

class TSGraph(object):
    '''
    Wrapper to represent a database object as a graph
    
    This class is primarily used by DoubleEndedConnect and has a number
    of utility functions for that purpose.  If you want access to the
    networkx Graph object, it is stored as self.graph
    
    Parameters
    ----------
    database :
        the database object to represent
    minima : list of minima, optional
        if none, include all minima and transition states from the database, 
        else include only the minima in the list
    no_edges : bool, optional
        If true include no transition states in the graph.
        This is used, for example, to try to find a new connection between 
        two minima.
    
    See Also
    --------
    DoubleEndedConnect
    
    Examples
    --------
    a graph can be easily constructed from a database::
    
    >>> graph = TSGraph(database)
    
    the networkx graph is accessed directly by
    
    >>> networkx_graph = graph.graph
    '''
    
    def __init__(self, database, minima=None, no_edges=False):
        self.graph = nx.Graph()
        self.storage = database
        self.connected_components = _ConnectedComponents(self.graph)
        self.minima = minima
        self.no_edges = no_edges
        self.refresh()
    
    def _build_all(self):
        """
        add all minima and all transition states to the graph
        """
        for m in self.storage.minima():
            self.graph.add_node(m)
        if not self.no_edges:
            for ts in self.storage.transition_states():
                self.graph.add_edge(ts.minimum1, ts.minimum2, ts=ts)

    def _build_from_list(self, minima):
        """
        add only those minima from the list `minima` to the graph.
        Don't add any transition states
        """
        minima = set(minima)
        for m in minima:
            self.graph.add_node(m)
        if not self.no_edges:       
            for ts in self.storage.transition_states():
                m1, m2 = ts.minimum1, ts.minimum2
                if m1 in minima:
                    if m2 in minima:
                        self.graph.add_edge(ts.minimum1, ts.minimum2, ts=ts)

    def refresh(self):
        self.connected_components.setRebuild()
        if self.minima is None:
            self._build_all()
        else:
            self._build_from_list(self.minima)

    def addMinimum(self, energy, coords, **kwargs):
        """
        add a minimum to the database and graph
        """
        minimum = self.storage.addMinimum(energy, coords, **kwargs)
        if not self.graph.has_node(minimum):
            self.connected_components.setRebuild()
        self.graph.add_node(minimum)
        return minimum
    
    def addTransitionState(self, E, coords, min1, min2, **kwargs):
        self.connected_components.setRebuild()
        ts = self.storage.addTransitionState(E, coords, min1, min2, **kwargs)
        self.graph.add_edge(ts.minimum1, ts.minimum2, ts=ts)
        return ts
            
    def areConnected(self, min1, min2):
        return self.connected_components.areConnected(min1, min2)
        try:
            ret = nx.bidirectional_dijkstra(self.graph, min1, min2)
            if ret != None:
                return True
            else: 
                return False
        except nx.NetworkXNoPath:
            return False
        
    def getPath(self, min1, min2):
        try:
            return nx.bidirectional_dijkstra(self.graph, min1, min2)
        except nx.NetworkXNoPath:
            return None

    def mergeMinima(self, min1, min2, update_database=True):
        """
        delete minima2.  all transition states pointing to min2 should
        now point to min1
        """
        self.connected_components.setRebuild()
        #make the edges of min2 now point to min1
        for v, data in self.graph[min2].iteritems():
            if v == min1: continue
            if v == min2: 
                print "warning: minimum", min2._id, "is connected to itself"
                continue
            #the new edge will be (min1, v).  Add it if it doesn't already exist
            if not self.graph.has_edge(min1, v):
#                if not self.graph.has_edge(v, min1):
#                    data = self.graph.get_edge_data(min2, v)
                self.graph.add_edge(min1, v, **data)
        self.graph.remove_node(min2)

        if update_database:
            self.storage.mergeMinima(min1, min2)



#
# below here only for testing
#

def create_random_database(nmin=20, nts=None, natoms=2):
    """
    create a database for test purposes
    """
    from pygmin.storage import Database
    import numpy as np
    
    if nts is None:
        nts = nmin
    db = Database()
    #generate random structures
    minlist = []
    for i in range(nmin):
        coords = np.random.uniform(-1,1,natoms*3)
        e = float(i) #make up a fake energy
        minlist.append( db.addMinimum(e, coords) )
    #add random transition states
    for i in range(nts):
        j1, j2 = 1, 1
        while j1 == j2:
            j1, j2 = np.random.randint(0, nmin, 2)
        m1, m2 = minlist[j1], minlist[j2] 
        coords = np.random.uniform(-1,1,natoms*3)
        e = float(j1 + j2)
        db.addTransitionState(e, coords, m1, m2)
    return db

class Graph(TSGraph):
    """this is included for backwards compatibility"""
    pass



import unittest
class TestGraph(unittest.TestCase):
    def setUp(self):
        self.db = create_random_database()
    
    def test_hi(self):
        graph = TSGraph(self.db)
        ng = graph.graph.number_of_nodes()
        nd = len(self.db.minima()) 
        self.assertEqual(ng, nd, "all nodes not imported")
        
        ntsg = graph.graph.number_of_edges()
        ntsd = len(self.db.transition_states()) 
        self.assertEqual(ntsg, ntsd, "all transition states not imported")
    
    def test_minima_keyword(self):
        nmin = 3
        minima = list(self.db.minima())
        minima = minima[:nmin]
        graph = TSGraph(self.db, minima=minima, no_edges=True)
        ng = graph.graph.number_of_nodes()
        self.assertEqual(ng, nmin)
        
        self.assertEqual(graph.graph.number_of_edges(), 0)
    
    def test_merge_minima(self):
        graph = TSGraph(self.db)
        #select two minima with at least one edge
        minima = []
        for n in graph.graph.nodes():
            if graph.graph.degree(n) > 0:
                minima.append(n)
            if len(minima) == 2:
                break
        nedges = [graph.graph.degree(n) for n in minima]
        min1, min2 = minima
        neighbors = graph.graph.neighbors(min2)
        graph.mergeMinima(min1, min2, update_database=True)
        self.assertNotIn(min2, graph.graph.nodes())
        
        #make sure the edges were copied
        for n in neighbors:
            if n == min1: continue
            edge_exists = (graph.graph.has_edge(min1, n)
                           or graph.graph.has_edge(n, min1))
            self.assert_(edge_exists)
            

        #make sure the new edges have transition state data attached
        for v in graph.graph.neighbors(min1):
            data = graph.graph.get_edge_data(v, min1)
            self.assertIn("ts", data.keys())

        #make sure min2 is not in database
        self.assertNotIn(min2, self.db.minima())
        
    def test_networkx(self):
        """check how networkx works"""
        graph = nx.Graph()
        graph.add_edge(1, 2, ts=1)
        graph.add_edge(2, 1, ts=2)
        self.assertEqual(graph.number_of_edges(), 1)
        
        data = graph.get_edge_data(1, 2)
        self.assertIsNotNone(data)
        self.assertDictEqual(graph.get_edge_data(1, 2),
                             graph.get_edge_data(2, 1))
        
        self.assertTrue(graph.has_edge(1,2))
        self.assertTrue(graph.has_edge(2,1))

if __name__ == "__main__":
    unittest.main()

