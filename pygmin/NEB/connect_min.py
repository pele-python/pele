import numpy as np
from pygmin.NEB.NEB import NEB
from pygmin.optimize.transition_state.transition_state_refinement import findTransitionState
import tstools
from pygmin.storage.savenlowest import Minimum 
from min_ts_graph import TransitionState
import networkx as nx
import copy
import itertools


class DoubleEndedConnect(object):
    def __init__(self, min1, min2, pot, graph, mindist):
        self.graph = graph
        self.minstart = min1
        self.minend = min2
        self.minstartid = graph.addMin(self.minstart)
        self.minendid = graph.addMin(self.minend)
        self.pot = pot
        self.mindist = mindist
        self.idlist = []
        self.distmatrix = dict()
        self.pairsNEB = dict()
        
    def _getDist(self, min1id, min2id):
        min1 = self.graph.node2min[min1id]
        min2 = self.graph.node2min[min2id]
        return self.getDist(min1, min2)
    def getDist(self, min1, min2):
        id1 = self.graph.minToKey(min1)
        id2 = self.graph.minToKey(min2)
        dist = self.distmatrix.get((id1, id2))
        if dist == None:
            dist = self.distmatrix.get((id2, id1))
        if dist != None:
            return dist
        print "calculating distance between", id1, id2
        dist, coords1, coords2 = self.mindist(min1.coords, min2.coords)
        self.distmatrix[(id1,id2)] = dist
        return dist

    def addTS(self, min1, min2, ts):
        id1 = self.addMin(min1)
        id2 = self.addMin(min2)
        if id1 == id2:
            print "warning: stepping off the transition state resulted in twice the same minima", id1
        else:
            print "adding transition state", id1, id2
            self.graph.addTS(min1, min2, ts)

    def addMin(self, min):
        id = self.graph.addMin(min)
        if not id in self.idlist and id != self.minstartid and id != self.minendid:
            self.idlist.append(id)
        return id

    def doNEB(self, min1, min2):
        graph = self.graph
        min1id = graph.minToKey(min1)
        min2id = graph.minToKey(min2)
        self.pairsNEB[(min1id, min2id)] = True
        print ""
        print "starting NEB run to try to connect minima", min1id, min2id
        print "minimizing the distance between the two minima"
        dist, newcoords1, newcoords2 = self.mindist(min1.coords, min2.coords) 
        print "runing NEB"
        neb = NEB(newcoords1, newcoords2, self.pot)
        neb.optimize()
        neb.MakeAllMaximaClimbing()
        neb.optimize()
    
        nclimbing = np.sum( [neb.isclimbing[i] for i in range(neb.nimages) ])
        print "from NEB search found", nclimbing, "climbing images"
        
        for i in range(neb.nimages):
            if neb.isclimbing[i]:
                coords = neb.coords[i,:]
                print "refining transition state from NEB climbing image"
                ret = findTransitionState(coords, self.pot)
                coords, eigval, eigvec, E, grad, rms = ret
    
                ts = TransitionState(coords, E, eigvec, eigval)
                print "falling off either side of transition state to find new minima"
                ret1, ret2 = tstools.minima_from_ts(self.pot.getEnergyGradient, coords, n = ts.eigenvec, \
                    displace=1e-3, quenchParameters={"tol":1e-7})
                
                newmin1 = Minimum(ret1[1], ret1[0])
                newmin2 = Minimum(ret2[1], ret2[0])
                self.addTS( newmin1, newmin2, ts)
                
                connected = graph._areConnected(min1id, min2id)
                print "connected yet?", connected
                if connected:
                    break
                
    def getNextPair(self):
        """
        this is the function which attempts to find a clever pair of minima to try to connect
        with the ultimate goal of connecting minstart and minend
        """
    
        #sort minimalist by distance to min1
        print "getting next pair to try to connect"
        print "sorting minimalist"
        minimalistbackup = copy.copy(self.idlist)
        def getDist(id1):
            min = self.graph.node2min[id1]
            return self.getDist(self.minstart, min)
        self.idlist = sorted(self.idlist, key=getDist)
    
        idlist = copy.copy(self.idlist)
        idlist.insert(0, self.minstartid)
        idlist.append(self.minendid)
        print idlist
        
        def mydist( mpair ):
            return self._getDist(mpair[0], mpair[1])
        pairlist = sorted( itertools.combinations(idlist,2), key=mydist )

        if True:
            
            if True:
                for id1, id2 in pairlist:
                    print id1, id2, self._getDist( id1, id2 )

            for id1, id2 in pairlist:
                if self.graph._areConnected(id1, id2):
                    continue
                if self.pairsNEB.has_key( (id1,id2) ) or self.pairsNEB.has_key( (id2,id1) ):
                    continue
                min1 = self.graph.node2min[id1]
                min2 = self.graph.node2min[id2]
                return min1, min2

        if False:
            #we should choose the next minima pair to try to connect in a more clever way.
            for i in range(len(idlist)-1):
                id1 = idlist[i]
                id2 = idlist[i+1]
                if self.graph._areConnected(id1, id2):
                    continue
                if self.pairsNEB.has_key( (id1,id2) ) or self.pairsNEB.has_key( (id2,id1) ):
                    continue
                min1 = self.graph.node2min[id1]
                min2 = self.graph.node2min[id2]
                return min1, min2
    
    
        return None, None
    
    def connect(self):
        """
        the main loop of the algorithm
        """
        self.doNEB(self.minstart, self.minend)
        while True: 
            if self.graph._areConnected(self.minstartid, self.minendid):
                return
            
            min1, min2 = self.getNextPair()
            if min1 == None:
                break
            self.doNEB(min1, min2)
        print "failed to find connection between", self.minstartid, self.minendid
    
    


def test():
    from min_ts_graph import MinTSGraph, getSetOfMinLJ
    from pygmin.optimize.quench import lbfgs_py as quench
    from pygmin.mindist.minpermdist_stochastic import minPermDistStochastic as mindist
    natoms = 27
    #get min1
    pot, saveit = getSetOfMinLJ(natoms)
    min1 = saveit.data[0]
    min2 = saveit.data[1]
    
    graph = MinTSGraph()
    min1id = graph.addMin(min1)
    min2id = graph.addMin(min2)
 
    connect = DoubleEndedConnect(min1, min2, pot, graph, mindist)
    connect.connect()
    
    print graph
    for node in graph.graph.nodes():
        print node, graph.node2min[node].E
    for edge, ts in graph.tstates().items():
        node1, node2 = edge 
        E1 = graph.node2min[node1].E
        E2 = graph.node2min[node2].E
        Ets = ts.E
        print node1, node2, "E", E1, Ets, E2
        
    ret = graph.getPath(min1, min2)
    if ret == None:
        print "no path found"
        return
    distances, path = ret
    with open("path.out", "w") as fout:
        for i in range(len(path)-1):
            n1 = path[i]
            n2 = path[i+1]
            m1 = graph.node2min[n1]
            m2 = graph.node2min[n2]
            ts = graph._getTS(n1, n2)
            print "path", n1, "->", n2, m1.E, "/->", ts.E, "\->", m2.E
            fout.write("%f\n" % m1.E)
            fout.write("%f\n" % ts.E)
        n2 = path[-1]
        m2 = graph.node2min[n2]
        fout.write("%f\n" % m2.E)


if __name__ == "__main__":
    test()


    