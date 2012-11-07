'''Wrapper to represent a storage class as a graph'''
import networkx as nx

__all__ = ["Graph"]


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
            

class Graph(object):
    '''
    Wrapper to represent a database object as a graph
    
    Parameter
    ---------
    database :
        the database object to represent
    
    Examples
    --------
    >>> graph = Graph(mydatabase)
    '''
    
    def __init__(self, database):
        self.graph=nx.Graph()
        self.storage = database
        self.connected_components = _ConnectedComponents(self.graph)
        self.refresh()
        
    def refresh(self):
        self.connected_components.setRebuild()
        for m in self.storage.minima():
            self.graph.add_node(m)
        for ts in self.storage.transition_states():
            self.graph.add_edge(ts.minimum1, ts.minimum2, ts=ts)

    def addMinimum(self, min1, *args):
        if not self.graph.has_node(min1):
            self.connected_components.setRebuild()
        return self.storage.addMinimum(min1, *args)
    
    def addTransitionState(self, *args, **kwargs):
        self.connected_components.setRebuild()
        return self.storage.addTransitionState(*args, **kwargs)
            
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
