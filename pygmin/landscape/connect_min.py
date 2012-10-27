import numpy as np
import networkx as nx
import copy
import itertools

from pygmin.transition_states import NEB, InterpolatedPath, findTransitionState, minima_from_ts
import pygmin.defaults as defaults
from pygmin.landscape import Graph

__all__ = ["DoubleEndedConnect"]

class DoubleEndedConnect(object):
    """
    Find a connected network of minima and transition states between min1 and min2
    
    Parameters
    ----------
    min1, min2 : Mimumum() objects
        the two minima to try to connect
    pot : potential object
        the potential
    mindist : callable
        the function which returns the optimized minimum distance between
        two structures
    database : Database() object
        the database object, used to save distance calculations so
        mindist() need only be called once for each minima pair. *Note* the
        use of this and graph is a bit redundant, this should be cleaned up
    tsSearchParams: dict
        parameters passed to the transition state search algorithm
    NEBquenchParams : dict
        parameters passed to the NEB minimization routine
    use_all_min : bool
        if True, then all known minima and transition states in graph will
        be used to try to connect min1 and min2.  This requires a mindist()
        call (or a retrieveal operation from database) for every pair which
        can take a very long time if many minima are known.  
    verbosity : int
        this controls how many status messages are printed.  (not really
        implemented yet)
    NEB_image_density : float
        how many NEB images per unit distance to use.
    NEB_iter_density : float
    NEBparams : dict
        NEB setup parameters.  E.g. this is used to pass the spring constant.
        (note: this is not for parameters related to interpolation).
    nrefine_max : int
        the maximum number of NEB transition state candidates to refine
    reoptimize_climbing : int
        the number of iterations to use for re-optimizing the climbing images
        after the NEB is done.
    
    Notes
    -----
    The algorithm is iterative, with each iteration composed of
    
    While min1 and min2 are not connected:
        1) choose a pair of know minima to try to connect
        
        2a) use NEB to get a guess for the transition state between them
        
        2b) refine the transition state to desired accuracy
        
        3) fall off either side of the transition state to find the two
        minima associated with that candidate
        
        4) add the transition state and associated minima to the known
        network
        
    
    Of the above, steps 1 and 2 are the most involved.  See the NEB and
    FindTransitionState classes for detailed descriptions of those
    routines.
    
    An important note is that the NEB is used only to get a *guess* for the
    transition state.  Thus we only want to put enough time and energy into
    the NEB routine to get the guess close enough that FindTransitionState
    can refine it to the correct transition state.  FindTransitionState is
    very fast if the initial guess is good, but can be very slow otherwise.
    
    Choose a pair
    -------------
    Here I will describe step 1), the algorithm to find a pair of known
    minima to try to connect.  This choice will keep in mind that the
    ultimate goal is to connect min1 and min2.
    
    In addition to the input parameter "graph", we keep a second graph
    "Gdist" which also has minima as the vertices. Gdist has an edge
    between (almost) every minima. The edge weight between vertices u and v
    is
    
        if u and v are connected in "graph":
            weight(u,v) = 0.
        else:
            weight(u,v) = mindist(u,v)
    
    Also, edges are removed from Gdist when an NEB is done to try to
    connect them.  This is to ensure we don't repeat NEB runs over and over
    again.  The minimum weight path between min1 and min2 in Gdist gives a
    good guess for the best way to try connect min1 and min2.  So the
    algorithm to find a pair of know minima (trial1, trial2) to try to
    connect is 
    
    path = Gdist.minimum_weight_path(min1, min2)
    trial1, trial2 = minima pair in path with lowest nonzero edge weight

    todo:
        use_all_min option is broken.  Fix up graph creation in general
        
        allow user to pass graph
        
        fail gracefully
        
    """    
    def __init__(
             self, min1, min2, pot, mindist, database, tsSearchParams=dict(), 
             NEBquenchParams = dict(), use_all_min=False, verbosity=1,
             NEB_image_density = 10., NEB_iter_density=15., NEBparams=dict(), 
             nrefine_max=100, reoptimize_climbing=0):
        self.minstart = min1
        self.minend = min2
        self.pot = pot
        self.mindist = mindist
        self.distmatrix = dict()
        self.pairsNEB = dict()
        self.idlist = []
        self.tsSearchParams = tsSearchParams
        self.NEBquenchParams = NEBquenchParams
        self.database = database
        self.graph = Graph(self.database)
        self.verbosity = int(verbosity)
        self.nrefine_max = nrefine_max
        
        self.NEB_image_density = float(NEB_image_density)
        self.NEB_iter_density = float(NEB_iter_density)
        self.NEBparams = NEBparams
        self.reoptimize_climbing = reoptimize_climbing

        self.Gdist = nx.Graph() 
        self._initializeDistances()
        self._initializeGdist(use_all_min)
        
        print "************************************************************"
        print "staring a double ended connect run between"
        print "        minimum 1: id %d energy %f" % (self.minstart._id, self.minstart.energy)
        print "        minimum 2: id %d energy %f" % (self.minend._id, self.minend.energy)
        print "        dist %f" % self.getDist(self.minstart, self.minend)
        print "************************************************************"
        
    
    def _initializeGdist(self, use_all_min):
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
                self._addNodeGdist(m)

    def _initializeDistances(self):
        """put all distances in the database into distmatrix for faster access"""
        for d in self.database.distances():
            id1 = d.minimum1._id
            id2 = d.minimum2._id
            
            if id2 > id1:
                id2, id1 = id1, id2
            self.distmatrix[(id1, id2)] = d.dist
            #self.distmatrix[(id2, id1)] = d.dist


    def _addNodeGdist(self, min1):
        """
        add node to graph Gdist.  Add an edge to all existing nodes with weight
        given by the distance.  If the minima are connected in self.graph make the 
        weight 0.
        """
        self.Gdist.add_node(min1)
        print "calculating distances between minimum", min1._id
        print "    and all others.  This could take some time"
        for min2 in self.Gdist.nodes():
            if min2 != min1:
                if not self.graph.areConnected(min1, min2):
                    dist = self.getDist(min1, min2)
                else:
                    dist = 0.
                self.Gdist.add_edge(min1, min2, {"dist":dist, "dist2":dist**2})

    def _remove_edgeGdist(self, min1, min2):
        try:
            self.Gdist.remove_edge(min1, min2)
        except nx.NetworkXError:
            pass
        return True


    def getDist(self, min1, min2):
        """
        get the distance between min1 and min2.

        get distances from the database if they exist, else calculate the
        distance and save it to the database
        """
        if True:
            id1, id2 = min1._id, min2._id
            if id2 > id1:
                id2, id1 = id1, id2
            dist = self.distmatrix.get((id1, id2))
            if dist is not None:
                return dist
        if self.database is None:
            return self._getDistNoDB(min1, min2)
        dist = self.database.getDistance(min1, min2)
        if dist is not None:
            return dist
        dist, coords1, coords2 = self.mindist(min1.coords, min2.coords)
        if self.verbosity > 1:
            print "calculated distance between", min1._id, min2._id, dist
        self.database.setDistance(dist, min1, min2)
        if True:
            id1, id2 = min1._id, min2._id
            if id2 > id1:
                id2, id1 = id1, id2
                self.distmatrix[(id1, id2)] = dist
        return dist

    
    def _getDistNoDB(self, min1, min2):
        dist = self.distmatrix.get((min1, min2))
        if dist is not None:
            return dist
        print "calculating distance between", min1._id, min2._id
        dist, coords1, coords2 = self.mindist(min1.coords, min2.coords)
        self.distmatrix[(min1,min2)] = dist
        self.distmatrix[(min2,min1)] = dist
        return dist

    def _addMinimum(self, *args):
        m = self.graph.addMinimum(*args)
        self.graph.refresh() #this could be slow
        if not m in self.idlist and not (m == self.minstart) and not (m == self.minend):
            self.idlist.append(m)
        self._addNodeGdist(m)
        return m    

    def _refineTS(self, coords):
        """
        find nearest transition state to NEB climbing image
        
        return True if we find a valid transition state
        
        return False if something goes wrong
        """
        #run ts search algorithm
        kwargs = dict(defaults.tsSearchParams.items() + self.tsSearchParams.items())
        ret = findTransitionState(coords, self.pot, **kwargs)
        
        #check to make sure it is a valid transition state 
        coords = ret.coords
        E = ret.energy
        if not ret.success:
            print "transition state search failed"
            return False
            
        if ret.eigenval >= 0.:
            print "warning: transition state has positive lowest eigenvalue, skipping:", ret.eigenval, ret.energy, ret.rms
            print "         not adding transition state"
            return False
        
        #find the minima which this transition state connects
        print "falling off either side of transition state to find new minima"
        ret1, ret2 = minima_from_ts(self.pot.getEnergyGradient, coords, n = ret.eigenvec, \
            displace=1e-3, quenchParameters={"tol":1e-7})
        
        #add the new minima to the graph (they may already be in there)
        min1 = self._addMinimum(ret1[1], ret1[0])
        min2 = self._addMinimum(ret2[1], ret2[0])
        if min1 == min2:
            print "warning: stepping off the transition state resulted in twice the same minima", min1._id
            return False
        
        #the transition state is good, add it to the graph
        print "adding transition state", min1._id, min2._id
        ts = self.graph.addTransitionState(E, coords, min1, min2, eigenvec=ret.eigenvec, eigenval=ret.eigenval)
        self.graph.refresh()
        """
        set the weight to zero in Gdist.
        """
        self.Gdist.add_edge(min1, min2, {"dist":0., "dist2":0.})
        return True
        
   
    def _doNEB(self, minNEB1, minNEB2, repetition = 0):
        """
        do NEB between minNEB1 and minNEB2.  refine any transition state candidates and
        minimize from either side of the transition state to find two new minima.
        """
        graph = self.graph
        
        #Make sure we haven't already tried this NEB pair and
        #record some data so we don't try it again in the future
        if self.pairsNEB.has_key((minNEB1, minNEB2)) and repetition == 0:
            print "WARNING: redoing NEB for minima", minNEB1._id, minNEB2._id
            print "         aborting NEB"
            self._remove_edgeGdist(minNEB1, minNEB2)
            return True
        self.pairsNEB[(minNEB1, minNEB2)] = True
        self.pairsNEB[(minNEB2, minNEB1)] = True

        #arrange the coordinates to minimize the distance between them        
        dist, newcoords1, newcoords2 = self.mindist(minNEB1.coords, minNEB2.coords)
        print ""
        
        if repetition == 0: 
            factor = 1.
        else: 
            #change parameters for second repetition
            print "running NEB a second time"
            print "    doubling the number of images"
            print "    doubling the number of steps"
            factor = float(repetition + 1)
        
        #determine the number of images
        nimages = int(max(1., dist) * self.NEB_image_density * factor)
        
        #determine the number of iterations
        NEBquenchParams = copy.copy(self.NEBquenchParams)
        if NEBquenchParams.has_key("nsteps"):
            niter = NEBquenchParams["nsteps"]
        else:
            niter = int(self.NEB_iter_density * nimages)
            NEBquenchParams["nsteps"] = niter
        
        
        #run NEB 
        NEBparams = dict(defaults.NEBparams.items() + self.NEBparams.items())
        print "starting NEB run to try to connect minima", minNEB1._id, minNEB2._id, dist
        print "    nimages", nimages
        print "    nsteps ", niter
        neb = NEB(InterpolatedPath(newcoords1, newcoords2, nimages), 
                  self.pot, **NEBparams)
        neb.optimize(**NEBquenchParams)
        neb.MakeAllMaximaClimbing()

        if self.reoptimize_climbing > 0:
            print "optimizing climbing images for a small number of steps"
            NEBquenchParams["nsteps"] = self.reoptimize_climbing
            neb.optimize(**NEBquenchParams)

    
        #get the transition state candidates from the NEB result
        climbing_images = [ (neb.energies[i], i) for i in range(neb.nimages) 
                           if neb.isclimbing[i] ]
        nclimbing = len(climbing_images)
        print "from NEB search found", nclimbing, "transition state candidates"
        if nclimbing == 0:
            print "WARNING: found zero climbing images.  Are the minima really the same?"
            print "         energies:", minNEB1.energy, minNEB2.energy, "distance", dist  
            return False   
        climbing_images = sorted(climbing_images, reverse=True) #highest energies first

        #find the nearest transition state the the transition state candidate
        success = False
        
        nrefine = min(self.nrefine_max, len(climbing_images))
        count = 0
        for energy, i in climbing_images[:nrefine]:
            count += 1
            print "refining transition state from NEB climbing image:", count, "out of", nrefine
            coords = neb.coords[i,:]
            ts_success = self._refineTS(coords)
            if ts_success:
                success = True
        
        #if the minima are still not connected, remove this edge so we don't try this NEB again
        if not self.graph.areConnected(minNEB1, minNEB2):
            self._remove_edgeGdist(minNEB1, minNEB2)
            
        return success
            
                
    def _getNextPair(self):
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
        #weight the edge by the square of the distance
        #could also be "dist"
        pathweight = "dist2" 
        try:
            path = nx.shortest_path(
                    self.Gdist, self.minstart, self.minend, weight=pathweight)
        except nx.NetworkXNoPath:
            print "Can't find any way to try to connect the minima"
            return None, None
        
        print "best guess for path.  (dist=0.0 means the path is known)"
        
        #path is a list of nodes
        distances = nx.get_edge_attributes(self.Gdist, pathweight)
        distmin = 1e100
        minpair = (None, None)
        for i in range(1,len(path)):
            min1 = path[i-1]
            min2 = path[i]
            dist = distances.get((min1,min2))
            if dist is None:
                dist = distances.get((min2,min1))
            #print "dist", dist
            if pathweight == "dist2":
                dist = np.sqrt(dist)
            print "    path guess", min1._id, min2._id, dist
            if dist > 1e-10 and dist < distmin:
                distmin = dist
                minpair = (min1, min2)
        return minpair
                
    
    def connect(self, maxiter=100):
        """
        the main loop of the algorithm
        """
        NEBattempts = 2;
        for i in range(maxiter): 
            if self.graph.areConnected(self.minstart, self.minend):
                return
            
            print ""
            print "starting connect cycle", i
            min1, min2 = self._getNextPair()
            if min1 is None or min2 is None:
                break
            for i in range(NEBattempts):
                NEB_success = self._doNEB(min1, min2, i)
                if NEB_success:
                    break
        print "failed to find connection between", self.minstart._id, self.minend._id

    def returnPath(self):
        """return information about the path"""
        if not self.graph.areConnected(self.minstart, self.minend):
            return
        minima = nx.shortest_path(self.graph.graph, self.minstart, self.minend)
        transition_states = []
        mints = [minima[0]]
        for i in range(1,len(minima)):
            m1 = minima[i-1]
            m2 = minima[i]
            ts = self.database.getTransitionState(m1, m2)
            transition_states.append(ts)
            mints.append(ts)
            mints.append(m2)
        
        S = np.zeros(len(mints))
        for i in range(1,len(mints)):
            coords1 = mints[i-1].coords
            coords2 = mints[i].coords
            dist, c1, c2 = self.mindist(coords1, coords2)
            S[i] = S[i-1] + dist
        energies = np.array([m.energy for m in mints])
        return mints, S, energies
        


###########################################################
#only testing stuff below here
###########################################################

def getSetOfMinLJ(natoms = 11): #for testing purposes
    from pygmin.potentials.lj import LJ
    pot = LJ()
    coords = np.random.uniform(-1,1,natoms*3)
    from pygmin.basinhopping import BasinHopping
    from pygmin.takestep.displace import RandomDisplacement
    from pygmin.takestep.adaptive import AdaptiveStepsize
    from pygmin.storage.database import Database
    import os
    #dbfile = "test.db"
    #os.remove(dbfile)
    #saveit = Database(db=dbfile)
    saveit = Database()
    takestep1 = RandomDisplacement()
    takestep = AdaptiveStepsize(takestep1, frequency=15)
    bh = BasinHopping(coords, pot, takestep, storage=saveit.minimum_adder(), outstream=None)
    bh.run(100)
    return pot, saveit


def test():
    from pygmin.landscape import Graph
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
    minima = database.minima()
    min1 = minima[0]
    min2 = minima[1]
    print min1.energy, min2.energy
    
    if False:
        #test to see if min1 and min2 are already connected
        connected = graph.areConnected(min1, min2)
        print "at start are minima connected?", connected
        return
 
    connect = DoubleEndedConnect(min1, min2, pot, mindist, database)
    connect.connect()
    
    graph = connect.graph
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


    
