import numpy as np
import networkx as nx

class TransitionState(object):
    """
    storage class for transition states
    """
    def __init__(self, coords, E, eigenvec, eigenval):
        self.coords = np.copy(coords)
        self.E = E
        self.eigenvec = np.copy(eigenvec)
        self.eigenval = eigenval


def sameMin(min1, min2, Ediff = 1e-3):
    return np.abs(min1.E - min2.E) <= Ediff

class MinTSGraph(object):
    """
    defines the graph object where the nodes represent minima and the edges represent 
    transition states.
    
    attached to each node is a Minimum object
    
    attached to each edge is a TransitionState object
    """
    def __init__(self):
        self.graph = nx.Graph()
        self.node2min = dict()
        self.edge2ts = dict()
        nx.set_edge_attributes(self.graph, "ts", dict())
        #nx.set_node_attributes( self.graph )
        self.nmin = 0
        #self.Etol = 1e-3
    
    def _newKey(self):
        self.nmin += 1
        return self.nmin

    def minToKey(self, min1):
        for node in self.graph.nodes():
            min2 = self.node2min[node]
            if sameMin(min1, min2):
                return node
        return None
    
    def addMin(self, min):
        key = self.minToKey(min)
        if key != None: return key
        
        #it's not in the graph, so get a label for it
        key = self._newKey()
        self.graph.add_node(key)
        self.node2min[key] = min
        return key
    
    def addTS(self, min1, min2, ts):
        minkey1 = self.minToKey(min1)
        minkey2 = self.minToKey(min2)
        return self._addTS(minkey1, minkey2, ts)
    
    def _addTS(self, minkey1, minkey2, ts):
        self.graph.add_edge(minkey1, minkey2, ts=ts)
    
    def areConnected(self, min1, min2):
        minkey1 = self.minToKey(min1)
        minkey2 = self.minToKey(min2)
        return self._areConnected(minkey1, minkey2)
    
    def _areConnected(self, minkey1, minkey2):
        try:
            ret = nx.bidirectional_dijkstra(self.graph, minkey1, minkey2)
            if ret != None:
                return True
            else: return False
        except nx.NetworkXNoPath:
            return False
        
    def getPath(self, min1, min2):
        minkey1 = self.minToKey(min1)
        minkey2 = self.minToKey(min2)
        return self._getPath(minkey1, minkey2)
    def _getPath(self, minkey1, minkey2):
        try:
            return nx.bidirectional_dijkstra(self.graph, minkey1, minkey2)
        except nx.NetworkXNoPath:
            return None

    
    def tstates(self):
        return nx.get_edge_attributes(self.graph, "ts")
    
    def getTS(self, min1, min2):
        minkey1 = self.minToKey(min1)
        minkey2 = self.minToKey(min2)
        return self._getTS(minkey1, minkey2)
        
    def _getTS(self, minkey1, minkey2):
        e2ts = self.tstates()
        try:
            return e2ts[(minkey1, minkey2)]
        except:
            return e2ts[(minkey2, minkey1)]

def getSetOfMinLJ(natoms = 11): #for testing purposes
    from pygmin.potentials.lj import LJ
    pot = LJ()
    coords = np.random.uniform(-1,1,natoms*3)
    from pygmin.basinhopping import BasinHopping
    from pygmin.takestep.displace import RandomDisplacement
    from pygmin.takestep.adaptive import AdaptiveStepsize
    from pygmin.storage.savenlowest import SaveN
    saveit = SaveN(10)
    takestep1 = RandomDisplacement()
    takestep = AdaptiveStepsize(takestep1, frequency=15)
    bh = BasinHopping(coords, pot, takestep, storage=saveit, outstream=None)
    bh.run(100)
    return pot, saveit

def test():
    from pygmin.optimize import lbfgs_py as quench
    from pygmin.storage.savenlowest import Minimum
    natoms = 11
    #get min1
    pot, saveit = getSetOfMinLJ(natoms)
    min1 = saveit.data[0]
    min2 = saveit.data[1]
    
    graph = MinTSGraph()
 
    node = graph.addMin(min1)
    print "node", node
    node = graph.addMin(min1)
    print "node", node
    node = graph.addMin(min2)
    print "node", node
    


if __name__ == "__main__":
    test()
