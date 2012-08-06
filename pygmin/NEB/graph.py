'''Wrapper to represent a storage class as a graph'''
import networkx as nx

class Graph(object):
    '''
    Wrapper to represent a storage class as a graph
    
    :example:
        
    >>> graph = Graph(save)
    '''
    storage=None
    graph=nx.Graph()
    
    def __init__(self, storage):
        '''
        :param storage: storage class to attach to        
        '''
        self.storage = storage
        self.refresh()
        
    def refresh(self):
        for m in self.storage.minima():
            self.graph.add_node(m)
        for ts in self.storage.transition_states():
            self.graph.add_edge(ts.minimum1, ts.minimum2, ts=ts)
            
    def addMinimum(self, *args):
        return self.storage.addMinimum(*args)
    
    def addTransitionState(self, *args, **kwargs):
        return self.storage.addTransitionState(*args, **kwargs)
            
    def areConnected(self, min1, min2):
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
        
        
    