import numpy as np
from pygmin.NEB.NEB import NEB
from pygmin.optimize.transition_state.transition_state_refinement import findTransitionState
import tstools
from pygmin.storage.savenlowest import Minimum 
import networkx as nx
import copy
import itertools


class DoubleEndedConnect(object):
    def __init__(
                 self, min1, min2, pot, graph, mindist, database=None, tsSearchParams = None, 
                 NEB_optimize_quenchParams = None, use_all_min=False):
        self.graph = graph
        self.minstart = min1
        self.minend = min2
        self.pot = pot
        self.mindist = mindist
        self.distmatrix = dict()
        self.pairsNEB = dict()
        self.idlist = []
        self.tsSearchParams = tsSearchParams
        self.NEB_optimize_quenchParams = NEB_optimize_quenchParams
        self.database = database
        

        self.Gdist = nx.Graph() 
        self.initializeGdist(use_all_min)
        
    
    def initializeGdist(self, use_all_min):
        """
        Initialize the graph Gdist. Gdist is used to help try to find the 
        optimal way to try to connect minstart and minend.
        
        Gdist has a vertex for minstart, minend and every new minima found.  (see note)
        
        Gdist has an edge between (almost) every minima.
        
        the edge weight is:
            if self.graph.areConnected(u,v):
                weight(u,v) = 0.
            else:
                weight(u,v) = self.getDist(u,v)
        
        edges are removed from the graph when an NEB is done to try to connect them.
        
        note: currently we initialize Gdist only with minstart and minend.  If we use
        all the other minima in self.graph we might be able to find a better path or
        find it more quickly.  But if the database of minima is large calculating the
        distances between every minima can take a long time
        """
        nx.set_edge_attributes(self.Gdist, "dist", dict())
        nx.set_edge_attributes(self.Gdist, "dist2", dict())
        dist = self.getDist(self.minstart, self.minend)
        self.Gdist.add_edge(self.minstart, self.minend, {"dist":dist, "dist2":dist**2})
        if use_all_min:
            """
            add all minima in self.graph to self.Gdist
            """
            for m in self.graph.graph.nodes():
                self.addNodeGdist(m)
    def addNodeGdist(self, min1):
        """
        add node to graph Gdist.  Add an edge to all existing nodes with weight
        given by the distance.  If the minima are connected in self.graph make the 
        weight 0.
        """
        self.Gdist.add_node(min1)
        for min2 in self.Gdist.nodes():
            if min2 != min1:
                if not self.graph.areConnected(min1, min2):
                    dist = self.getDist(min1, min2)
                else:
                    dist = 0.
                self.Gdist.add_edge(min1, min2, {"dist":dist, "dist2":dist**2})
                    

    def getDist(self, min1, min2):
        """
        get distances from the database if they exist, else calculate the
        distance and save it to the database
        """
        if self.database is None:
            return self.getDistNoDB(min1, min2)
        dist = self.database.getDistance(min1, min2)
        if dist is not None:
            return dist
        dist, coords1, coords2 = self.mindist(min1.coords, min2.coords)
        print "calculated distance between", min1._id, min2._id, dist
        self.database.setDistance(dist, min1, min2)
        #self.distmatrix[(min1,min2)] = dist
        #self.distmatrix[(min2,min1)] = dist
        return dist

    
    def getDistNoDB(self, min1, min2):
        dist = self.distmatrix.get((min1, min2))
        if dist is not None:
            return dist
        print "calculating distance between", min1._id, min2._id
        dist, coords1, coords2 = self.mindist(min1.coords, min2.coords)
        self.distmatrix[(min1,min2)] = dist
        self.distmatrix[(min2,min1)] = dist
        return dist

    def addMinimum(self, *args):
        m = self.graph.addMinimum(*args)
        self.graph.refresh() #this could be slow
        if not m in self.idlist and not (m == self.minstart) and not (m == self.minend):
            self.idlist.append(m)
        self.addNodeGdist(m)
        return m    
    
    def doNEB(self, minNEB1, minNEB2):
        """
        do NEB between minNEB1 and minNEB2.  refine any transition state candidates and
        minimize from either side of the transition state to find two new minima.
        """
        graph = self.graph
        #record some data so we don't try this NEB again
        if self.pairsNEB.has_key((minNEB1, minNEB2)):
            print "WARNING: redoing NEB for minima", minNEB1._id, minNEB2._id
        self.pairsNEB[(minNEB1, minNEB2)] = True
        self.pairsNEB[(minNEB2, minNEB1)] = True
        #self.Gdist.remove_edge(minNEB1, minNEB2)
        print ""
        #print "minimizing the distance between the two minima"
        dist, newcoords1, newcoords2 = self.mindist(minNEB1.coords, minNEB2.coords) 
        print "starting NEB run to try to connect minima", minNEB1._id, minNEB2._id, dist
        #print "runing NEB"
        neb = NEB(newcoords1, newcoords2, self.pot)
        neb.optimize(quenchParams = self.NEB_optimize_quenchParams)
        neb.MakeAllMaximaClimbing()
        #neb.optimize(quenchParams = self.NEB_optimize_quenchParams)
    
        nclimbing = np.sum( [neb.isclimbing[i] for i in range(neb.nimages) ])
        print "from NEB search found", nclimbing, "climbing images"
        if nclimbing == 0:
            print "WARNING: found zero climbing images.  Are the minima really the same?"
            print "         energies:", minNEB1.energy, minNEB2.energy, "distance", dist 
        
        climbing_images = [ (neb.energies[i], i) for i in range(neb.nimages) 
                           if neb.isclimbing[i] ]
        climbing_images = sorted(climbing_images, reverse=True) #highest energies first
        
        for energy, i in climbing_images:
            coords = neb.coords[i,:]
            np.savetxt("climbingimage", coords)
            #print "exiting", exit(1)
            print "refining transition state from NEB climbing image"
            ret = findTransitionState(coords, self.pot, tsSearchParams = self.tsSearchParams)
            #coords, eigval, eigvec, E, grad, rms = ret
            coords = ret.coords
            E = ret.energy
            if ret.eigenval >= 0.:
                print "warning: transition state has positive lowest eigenvalue, skipping:", ret.eigenval, ret.energy, ret.rms
            print "falling off either side of transition state to find new minima"
            ret1, ret2 = tstools.minima_from_ts(self.pot.getEnergyGradient, coords, n = ret.eigenvec, \
                displace=1e-3, quenchParameters={"tol":1e-7})
            
            min1 = self.addMinimum(ret1[1], ret1[0])
            min2 = self.addMinimum(ret2[1], ret2[0])
            if min1 == min2:
                print "warning: stepping off the transition state resulted in twice the same minima", min1._id
            else:
                print "adding transition state", min1._id, min2._id
                #ts = TransitionState(coords, E, eigvec, eigval)
                #self.graph.addTransitionState(min1, min2, ts)   self.addTS( newmin1, newmin2, ts)
                from pygmin.storage.database import TransitionState
                ts = self.graph.addTransitionState(E, coords, min1, min2, eigenvec=ret.eigenvec, eigenval=ret.eigenval)
                graph.refresh()
                """
                set the weight to zero in Gdist.
                """
#                    if (min1,min2) in self.Gdist.edges() or (min2,min1) in self.Gdist.edges():
#                    #path = nx.shortest_path(self.Gdist, min1, min2)
#                    #print "mypath", path
#                    #if path is not None and len(path) == 2:
#                        print "found edge"
                self.Gdist.add_edge(min1, min2, {"dist":0., "dist2":0.})
                
            connected = graph.areConnected(minNEB1, minNEB2)
            print "connected yet?", connected
            if connected:
                break
        #if the minima are still not connected, remove this edge so we don't try this NEB again
        if not self.graph.areConnected(minNEB1, minNEB2):
            self.Gdist.remove_edge(minNEB1, minNEB2)
            
                
    def getNextPair(self):
        """
        this is now obsolete
        
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
    
    def getNextPairNew(self):
        """
        this is the function which attempts to find a clever pair of minima to try to 
        connect with the ultimate goal of connecting minstart and minend
        
        this method can be described as follows:
        
        make a new graph Gnew which is complete (all vetices connected).  The edges
        have a weight given by
        
        if self.graph.areConnected(u,v):
            weight(u,v) = 0.
        else:
            weight(u,v) = mindist(u,v)
        
        if an NEB has been attempted between u and v then the edge is removed.
        
        we then find the shortest path between minstart and minend.  We return
        the pair in this path which has edge weight with the smallest non zero value 
        
        update: find the shortest path weighted by distance squared.  This penalizes finding
        the NEB between minima that are very far away.  (Does this too much favor long paths?)
        """
        print "finding a good pair to try to connect"
        try:
            path = nx.shortest_path(
                    self.Gdist, self.minstart, self.minend, weight="dist2")
        except nx.NetworkXNoPath:
            print "Can't find any way to try to connect the minima"
            return None, None
        
        print "best guess for path"
        
        #path is a list of nodes
        distances = nx.get_edge_attributes(self.Gdist, "dist")
        distmin = 1e100
        minpair = (None, None)
        for i in range(1,len(path)):
            min1 = path[i-1]
            min2 = path[i]
            dist = distances.get((min1,min2))
            if dist is None:
                dist = distances.get((min2,min1))
            #print "dist", dist
            print "    path guess", min1._id, min2._id, dist
            if dist > 1e-10 and dist < distmin:
                distmin = dist
                minpair = (min1, min2)
        return minpair
                
                
            
        
        

    
    def connect(self):
        """
        the main loop of the algorithm
        """
#        print "starting a double ended search to try to connect minima", self.minstart._id, self.minend._id
#        if self.graph.areConnected(self.minstart, self.minend):
#            print "aborting double ended connect:  minima are already connected"
#            return
#        
#        if True:
#            min1, min2 = self.getNextPairNew()
#
#        self.doNEB(self.minstart, self.minend)
        while True: 
            if self.graph.areConnected(self.minstart, self.minend):
                return
            
            min1, min2 = self.getNextPairNew()
            if min1 is None or min2 is None:
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
    from pygmin.storage.database import Database
    import os
    dbfile = "test.db"
    os.remove(dbfile)
    saveit = Database(db=dbfile)
    takestep1 = RandomDisplacement()
    takestep = AdaptiveStepsize(takestep1, frequency=15)
    bh = BasinHopping(coords, pot, takestep, storage=saveit.minimum_adder(), outstream=None)
    bh.run(100)
    return pot, saveit


def test():
    from graph import Graph
    from pygmin.optimize.quench import lbfgs_py as quench
    from pygmin.mindist.minpermdist_stochastic import minPermDistStochastic as mindist
    from pygmin.storage.database import Database
    import pygmin.defaults as defaults
    defaults.quenchParams = {"iprint": 1}
    natoms = 11
    #get min1
    pot, database = getSetOfMinLJ(natoms)
#    from pygmin.potentials.lj import LJ
#    pot = LJ()
#    saveit = Database(db="test.db")
    graph = Graph(database)
    minima = database.minima()
    min1 = minima[0]
    min2 = minima[1]
    print min1.energy, min2.energy
    
    if False:
        #test to see if min1 and min2 are already connected
        connected = graph.areConnected(min1, min2)
        print "at start are minima connected?", connected
        return
 
    connect = DoubleEndedConnect(min1, min2, pot, graph, mindist, database=database)
    connect.connect()
    
    if False:
        print graph
        for node in graph.graph.nodes():
            print node._id, node.energy
    for ts in graph.storage.transition_states():
        print ts.minimum1._id,ts.minimum2._id, "E", ts.minimum1.energy, ts.minimum2.energy, ts.energy
        
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


    