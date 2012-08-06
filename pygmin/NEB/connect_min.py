import numpy as np
from pygmin.NEB.NEB import NEB
from pygmin.optimize.transition_state.transition_state_refinement import findTransitionState
import tstools
from pygmin.storage.savenlowest import Minimum 
import networkx as nx
import copy
import itertools


class DoubleEndedConnect(object):
    def __init__(self, min1, min2, pot, graph, mindist):
        self.graph = graph
        self.minstart = min1
        self.minend = min2
        self.pot = pot
        self.mindist = mindist
        self.distmatrix = dict()
        self.pairsNEB = dict()
        self.idlist = []
        
    def getDist(self, min1, min2):
        dist = self.distmatrix.get((min1, min2))
        if dist == None:
            dist = self.distmatrix.get((min1, min2))
        if dist != None:
            return dist
        print "calculating distance between", min1._id, min2._id
        dist, coords1, coords2 = self.mindist(min1.coords, min2.coords)
        self.distmatrix[(min1,min2)] = dist
        return dist

    def addMinimum(self, *args):
        m = self.graph.addMinimum(*args)
        if not m in self.idlist and not (m is self.minstart) and not (m is self.minend):
            self.idlist.append(m)
        return m
        
    def doNEB(self, min1, min2):
        graph = self.graph
        self.pairsNEB[(min1, min2)] = True
        print ""
        print "starting NEB run to try to connect minima", min1._id, min2._id
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
                #coords, eigval, eigvec, E, grad, rms = ret
                coords = ret.coords
                E = ret.energy
                print "falling off either side of transition state to find new minima"
                ret1, ret2 = tstools.minima_from_ts(self.pot.getEnergyGradient, coords, n = ret.eigenvec, \
                    displace=1e-3, quenchParameters={"tol":1e-7})
                
                min1 = self.addMinimum(ret1[1], ret1[0])
                min2 = self.addMinimum(ret2[1], ret2[0])
                if min1 is min2:
                    print "warning: stepping off the transition state resulted in twice the same minima", min1._id
                else:
                    print "adding transition state", min1._id, min2._id
                    #ts = TransitionState(coords, E, eigvec, eigval)
                    #self.graph.addTransitionState(min1, min2, ts)   self.addTS( newmin1, newmin2, ts)
                    from pygmin.storage.database import TransitionState
                    ts = self.graph.addTransitionState(E, coords, min1, min2, eigenvec=ret.eigenvec, eigenval=ret.eigenval)
                    graph.refresh()
                    
                connected = graph.areConnected(min1, min2)
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
        def getDist(m):
            return self.getDist(self.minstart, m)
        self.idlist = sorted(self.idlist, key=getDist)
    
        idlist = copy.copy(self.idlist)
        idlist.insert(0, self.minstart)
        idlist.append(self.minend)
        print idlist
        
        def mydist( mpair ):
            return self.getDist(mpair[0], mpair[1])
        pairlist = sorted( itertools.combinations(idlist,2), key=mydist )

        if True:
            
            if True:
                for m1, m2 in pairlist:
                    print m1._id, m2._id, self.getDist( m1, m2 )

            for m1, m2 in pairlist:
                if self.graph.areConnected(m1, m2):
                    continue
                if self.pairsNEB.has_key( (m1,m2) ) or self.pairsNEB.has_key( (m2,m1) ):
                    continue
                return m1, m2

        if False:
            #we should choose the next minima pair to try to connect in a more clever way.
            for i in range(len(idlist)-1):
                m1 = idlist[i]
                m2 = idlist[i+1]
                if self.graph.areConnected(m1, m2):
                    continue
                if self.pairsNEB.has_key( (m1,m2) ) or self.pairsNEB.has_key( (m2,m1) ):
                    continue
                return m1, m2
    
    
        return None, None
    
    def connect(self):
        """
        the main loop of the algorithm
        """
        self.doNEB(self.minstart, self.minend)
        while True: 
            if self.graph.areConnected(self.minstart, self.minend):
                return
            
            min1, min2 = self.getNextPair()
            if min1 == None:
                break
            self.doNEB(min1, min2)
        print "failed to find connection between", self.minstart._id, self.minend._id
    
def getSetOfMinLJ(natoms = 11): #for testing purposes
    from pygmin.potentials.lj import LJ
    pot = LJ()
    coords = np.random.uniform(-1,1,natoms*3)
    from pygmin.basinhopping import BasinHopping
    from pygmin.takestep.displace import RandomDisplacement
    from pygmin.takestep.adaptive import AdaptiveStepsize
    from pygmin.storage.database import Storage
    saveit = Storage(db="test.db")
    takestep1 = RandomDisplacement()
    takestep = AdaptiveStepsize(takestep1, frequency=15)
    bh = BasinHopping(coords, pot, takestep, storage=saveit.minimum_adder(), outstream=None)
    bh.run(100)
    return pot, saveit


def test():
    from graph import Graph
    from pygmin.optimize.quench import lbfgs_py as quench
    from pygmin.mindist.minpermdist_stochastic import minPermDistStochastic as mindist
    from pygmin.storage.database import Storage
    natoms = 27
    #get min1
    pot, saveit = getSetOfMinLJ(natoms)
#    from pygmin.potentials.lj import LJ
#    pot = LJ()
#    saveit = Storage(db="test.db")
    graph = Graph(saveit)
    minima = saveit.minima()
    min1 = minima[0]
    min2 = minima[1]
    print min1.energy, min2.energy
 
    connect = DoubleEndedConnect(min1, min2, pot, graph, mindist)
    connect.connect()
    
    print graph
    for node in graph.graph.nodes():
        print node, node.energy
    for ts in graph.storage.transition_states():
        print ts.minimum1._id,ts.minimum2._id, "E", ts.minimum1.energy, ts.minimum2.energy, ts.minimum2.energy
        
    ret = graph.getPath(min1, min2)
    if ret == None:
        print "no path found"
        return
    distances, path = ret
    with open("path.out", "w") as fout:
        for i in range(len(path)-1):
            m1 = path[i]
            m2 = path[i+1]
            n1 = m1._id
            m2 = m2._id
#            ts = graph._getTS(n1, n2)
#            print "path", n1, "->", n2, m1.E, "/->", ts.E, "\->", m2.E
            fout.write("%f\n" % m1.energy)
            fout.write("%f\n" % ts.energy)
        m2 = path[-1]
        n2 = m2._id
        fout.write("%f\n" % m2.energy)


if __name__ == "__main__":
    test()


    