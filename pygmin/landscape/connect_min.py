import numpy as np
import networkx as nx
import copy
import itertools

from pygmin.transition_states import NEB, InterpolatedPath, findTransitionState, minima_from_ts
import pygmin.defaults as defaults
from pygmin.landscape import Graph, LocalConnect

__all__ = ["DoubleEndedConnect"]

class _DistanceGraph(object):
    """
    This graph is used to guess a good method for connecting two minima
    
    Parameters
    ----------
    database : 
        the database in which to store newly found minima and transition states.
        If database contains contains previously found minima and transition states,
        those will be used to help with the connection.
    graph :
        the graph build from the database which contains the minima and transition states
    mindist : callable
        the routine which calculates the optimized distance between two structures
    verbosity :
        how much info to print (not very thoroughly implemented)
    defer_database_update : bool
        if true, save new distances and only update the database when enough new
        distances have been accumulated
    db_update_min : int
        only update the database when at least this many new distances have been found.
    
    This graph has a vertex for every minimum and an edge
    between (almost) every minima pairs. The edge weight between vertices u and v
    is
    
        if u and v are connected in "graph":
            weight(u,v) = 0.
        else:
            weight(u,v) = mindist(u,v)
            (or  = mindist(u,v)**2)
    
    Also, edges are removed from this graph when an NEB is done to try to
    connect them.  This is to ensure we don't repeat NEB runs over and over
    again.  The minimum weight path between min1 and min2 in this graph gives a
    good guess for the best way to try connect min1 and min2.  So the
    algorithm to find a pair of know minima (trial1, trial2) to try to
    connect is 

    path = Gdist.minimum_weight_path(min1, min2)
    trial1, trial2 = minima pair in path with lowest nonzero edge weight

    """
    def __init__(self, database, graph, mindist, verbosity=0,
                 defer_database_update=True, db_update_min=300):
        self.database = database
        self.graph = graph
        self.mindist = mindist
        self.verbosity = verbosity
        
        self.Gdist = nx.Graph()
        self.distance_map = dict() #place to store distances locally for faster lookup
        nx.set_edge_attributes(self.Gdist, "weight", dict())
        self.debug = True
        
        self.defer_database_update = defer_database_update
        
        self.new_distances = dict() #keep track of newly calculated distances
        self.db_update_min = db_update_min

    def distToWeight(self, dist):
        """
        this function defines how the edge weight is calculated.
        
        good options might be:
        
        weight = dist
        
        weight = dist**2
            this favors paths with many short edges over
            paths with fewer but longer edges.
        """
        return dist**2

    def updateDatabase(self, force=False):
        """update databases with new distances"""
        nnewdist = len(self.new_distances.items())
        if nnewdist == 0:
            return
        if not force:
            if nnewdist < self.db_update_min:
                return
        print "updating database with", nnewdist, "new distances"
        self.database.setDistanceBulk(self.new_distances.iteritems())
        self.new_distances = dict()

    def _setDist(self, min1, min2, dist):
        """
        this function saves newly calculated distances both to the local
        distance map and ultimately to the database
        """
        #add the distance to the database and the distance map
        if self.defer_database_update:
            self.new_distances[(min1, min2)] = dist
        else:
            self.database.setDistance(dist, min1, min2)
        self.distance_map[(min1, min2)] = dist
        
        #make sure a zeroed edge weight is not overwritten
        #if not self.edge_weight.has_key((min1, min2)):
        #    if not self.edge_weight.has_key((min2, min1)):
        #        self.edge_weight[(min1, min2)] = weight
    
    def _getDistNoCalc(self, min1, min2):
        """
        get distance from local memory.  don't calculate it.
        """
        #first try to get the distance from the dictionary 
        dist = self.distance_map.get((min1,min2))
        if dist is not None: return dist
        dist = self.distance_map.get((min2,min1))
        if dist is not None: return dist

        if False:
            #this is extremely slow for large databases (> 50% of time spent here)
            #also, it's not necessary if we load all the distances in initialize()
            #if that fails, try to get it from the database
            dist = self.database.getDistance(min1, min2)
            if dist is not None: 
                print "distance in database but not in distance_map"
                return dist
        return None

    def getDist(self, min1, min2):
        """
        return the distance between two minima.  Calculate it and store it if
        not already known
        """
        dist = self._getDistNoCalc(min1, min2)
        if dist is not None: return dist
        
        #if it's not already known we must calculate it
        dist, coords1, coords2 = self.mindist(min1.coords, min2.coords)
        if self.verbosity > 1:
            print "calculated distance between", min1._id, min2._id, dist
        self._setDist(min1, min2, dist)
        return dist
    
    def _addEdge(self, min1, min2):
        """
        add a new edge to the graph.  Calculate the distance
        between the minima and set the edge weight
        """
        if min1 == min2: return
        dist = self.getDist(min1, min2)
        weight = self.distToWeight(dist)
        self.Gdist.add_edge(min1, min2, {"weight":weight})
        if self.graph.areConnected(min1, min2):
            self.setTransitionStateConnection(min1, min2)
            #note: this is incomplete.  if a new edge between
            #min1 and min2 connects
            #two clusters that were previously unconnected
            #then the edge weight should be set to zero 
            #with self.setTransitionStateConnection() for all minima
            #in the two clusters.  Currently this is being fixed
            #by calling checkGraph() from time to time.  I'm not sure
            #which is better.
    
    def _addMinimum(self, m):
        self.Gdist.add_node(m)
        #for noded that are connected set the edge weight using setTransitionStateConnection
        cc = nx.node_connected_component(self.graph.graph, m)
        for m2 in cc:
            if m2 in self.Gdist:
                #self.Gdist.add_edge(m, m2, weight=0.)
                self.setTransitionStateConnection(m, m2)
        
        #for all other nodes set the weight to be the distance
        for m2 in self.Gdist.nodes():
            if not self.Gdist.has_edge(m, m2):
                dist = self.getDist(m, m2)
                weight = self.distToWeight(dist)
                self.Gdist.add_edge(m, m2, {"weight":weight})


        
        
    
    def addMinimum(self, m):
        """
        add a new minima to the graph and add an edge to all the other
        minima in the graph.  
        
        Note: this can take a very long time if there are lots of minima
        in the graph.  mindist need to be run many many times.
        """
        trans = self.database.connection.begin()
        try:
            if not m in self.Gdist:
                self._addMinimum(m)
#                self.Gdist.add_node(m)
#                #add an edge to all other minima
#                for m2 in self.Gdist.nodes():
#                    self._addEdge(m, m2)
        except:
            trans.rollback()
            raise
        trans.commit()
                               

    def removeEdge(self, min1, min2):
        """remove an edge from the graph
        
        used to indicate that the routine should not try to connect
        these minima again.
        """
        self.Gdist.add_edge(min1, min2, weight=10e10)
        try:
            self.Gdist.remove_edge(min1, min2)
        except nx.NetworkXError:
            pass
        return True

    def _initializeDistances(self):
        """put all distances in the database into distance_map for faster access"""
#        from pygmin.storage.database import Distance
#        from sqlalchemy.sql import select
#        conn = self.database.engine.connect()
#        sql = select([Distance.__table__])
#        for tmp, dist, id1, id2 in conn.execute(sql):
#            #m1 = self.database.getMinimum(id1)
#            #m2 = self.database.getMinimum(id2)
#            self.distance_map[id1, id2] = dist
        if False:
            for d in self.database.distances():
                self.distance_map[(d.minimum1, d.minimum2)] = d.dist
        else:
            for d in self.database.distances():
                self.distance_map[(d._minimum1_id, d._minimum2_id)] = d.dist

    def replaceTransitionStateGraph(self, graph):
        self.graph = graph

    def _addRelevantMinima(self, minstart, minend):
        """
        add all the relevant minima from the database to the distance graph
        
        a minima is considered relevant if distance(min1, minstart) and
        distance(min1, minend) are both less than distance(minstart, minend)
        
        also, don't calculate any new distances, only add a minima if all distances
        are already known. 
        """
        start_end_distance = self.getDist(minstart, minend)
        count = 0
        naccept = 0
        for m in self.graph.graph.nodes():
            count += 1
            d1 = self._getDistNoCalc(m, minstart)
            if d1 is None: continue
            if d1 > start_end_distance: continue
            
            d2 = self._getDistNoCalc(m, minend)
            if d2 is None: continue
            if d2 > start_end_distance: continue
            
            print "    accepting minimum", d1, d2, start_end_distance
            
            naccept += 1
            self.addMinimum(m)
        print "    found", naccept, "relevant minima out of", count


    def initialize(self, minstart, minend, use_all_min=False, use_limited_min=True):
        """
        set up the distance graph
        
        initialize distance_map, add the start and end minima and load any other
        minima that should be used in the connect routine.
        """
        #raw_input("Press Enter to continue:")
        print "loading distances from database"
        self._initializeDistances()
        #raw_input("Press Enter to continue:")
        dist = self.getDist(minstart, minend)
        self.addMinimum(minstart)
        self.addMinimum(minend)
        if use_all_min:
            """
            add all minima in self.graph to self.Gdist
            """
            print "adding all minima to distance graph (Gdist)."
            print "    This might take a while."
            for m in self.database.minima():
                self.addMinimum(m)
        elif use_limited_min:
            print "adding relevant minima to distance graph (Gdist)."
            print "    This might take a while."
            self._addRelevantMinima(minstart, minend)
        #raw_input("Press Enter to continue:")

    def setTransitionStateConnection(self, min1, min2):
        """use this function to tell _DistanceGraph that
        there exists a known transition state connection between min1 and min2
        
        The edge weight will be set to zero
        """
        weight = 0.
        self.Gdist.add_edge(min1, min2, {"weight":weight})

    def shortestPath(self, min1, min2):
        """return the minimum weight path path between min1 and min2""" 
        if True:
            print "Gdist has", self.Gdist.number_of_nodes(), "nodes and", self.Gdist.number_of_edges(), "edges"
        try:
            path = nx.shortest_path(
                    self.Gdist, min1, min2, weight="weight")
        except nx.NetworkXNoPath:
            return None, None
        
        #get_edge attributes is really slow:
        #weights = nx.get_edge_attributes(self.Gdist, "weight") #this takes order number_of_edges
        weights = [ self.Gdist[path[i]][path[i+1]]["weight"] for i in range(len(path)-1) ]
        
        return path, weights
        
    def mergeMinima(self, min1, min2):
        """
        rebuild the graph with min2 deleted and 
        everything pointing to min1 pointing to min2 instead
        """
        print "    rebuilding Gdist"
        weights = nx.get_edge_attributes(self.Gdist, "weight")
        newgraph = nx.Graph()
        nx.set_edge_attributes(newgraph, "weight", dict())
        for node in self.Gdist.nodes():
            if node != min2:
                newgraph.add_node(node)
        for e in self.Gdist.edges():
            if not min1 in e and not min2 in e:
                newgraph.add_edge(e[0], e[1], {"weight":weights[e]})
            if min1 in e and min2 in e:
                continue
            #if e already exists in newgraph, make sure we don't overwrite
            #a zeroed edge weight
            if min2 in e:
                if e[0] == min2:
                    enew = (min1, e[1])
                else:
                    enew = (e[0], min1)
            else:
                enew = e
            existing_weight = weights.get(enew)
            if existing_weight is not None:
                if existing_weight < 1e-10:
                    #existing weight is zero.  don't overwrite
                    continue
            newgraph.add_edge(enew[0], enew[1], {"weight":weights[e]})
        
        #done, replace Gdist with newgraph
        self.Gdist = newgraph
        if self.debug:
            self.checkGraph()

    def checkGraph(self):
        """
        make sure graph is up to date.
        and make any corrections
        """
        print "checking Gdist"
        #check that all edges that are connected in self.graph
        #have zero edge weight
        #note: this could be done a lot more efficiently
        weights = nx.get_edge_attributes(self.Gdist, "weight")
        count = 0
        for e in self.Gdist.edges():
            are_connected = self.graph.areConnected(e[0], e[1])
            zero_weight = weights[e] < 1e-10
            #if they are connected they should have zero_weight
            if are_connected and not zero_weight:
                #print "    problem: are_connected", are_connected, "but weight", weights[e], "dist", dist
                if True:
                    #note: this is an inconsistency, but it's only a problem if
                    #there is no zero weight path from e[0] to e[1]
                    path, path_weight = self.shortestPath(e[0], e[1])
                    weight_sum = sum(path_weight)
#                    for i in range(len(path)-1):
#                        m1, m2 = path[i], path[i+1]
#                        try:
#                            w = weights[(m1, m2)]
#                        except KeyError:
#                            w = weights[(m2, m1)]
#                        weight_sum += w
                    if weight_sum > 10e-6:
                        #now there is definitely a problem.
                        count += 1
                        dist = self.getDist(e[0], e[1])
                        print "    problem: are_connected", are_connected, "but weight", weights[e], "dist", dist, "path_weight", weight_sum
                self.setTransitionStateConnection(e[0], e[1])
                            
                     
            if not are_connected and zero_weight:
                dist = self.getDist(e[0], e[1])
                print "    problem: are_connected", are_connected, "but weight", weights[e], "dist", dist
                w = self.distToWeight(dist)
                self.Gdist.add_edge(e[0], e[1], {"weight":w})
        if count > 0:
            print "    found", count, "inconsistencies in Gdist"
                     


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
    use_all_min : bool
        if True, then all known minima and transition states in graph will
        be used to try to connect min1 and min2.  This requires a mindist()
        call (or a retrieveal operation from database) for every pair which
        can take a very long time if many minima are known.  
    verbosity : int
        this controls how many status messages are printed.  (not really
        implemented yet)
    merge_minima : bool
        if True, minima for which NEB finds no transition state candidates 
        between them will be merged
    max_dist_merge : float
        merging minima will be aborted if the distance between them is greater
        than max_dist_merge
    local_connect_params : dict
        parameters passed to the local connect algorithm.  This includes all
        NEB and all transition state search parameters, along with, e.g. 
        now many times to retry a local connect run.  See documentation for
        LocalConnect for details.
    fresh_connect : bool
        if true, ignore all existing minima and transition states in the
        database and try to find a new path
    longest_first : bool
        if true, always try to connect the longest segment in the path guess
        first
    
    Notes
    -----
    The algorithm is iterative, with each iteration composed of
    
    While min1 and min2 are not connected:
        1) choose a pair of known minima to try to connect
        
        2) use NEB to get a guess for the transition states between them
        
        3) refine the transition states to desired accuracy
        
        4) fall off either side of the transition states to find the two
        minima associated with that candidate
        
        5) add the transition state and associated minima to the known
        network
        
    
    Of the above, steps 1 and 2 and 3 are the most involved.  See the NEB and
    FindTransitionState classes for detailed descriptions of steps 2 and 3.
    
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
    "Gdist" (now called "dist_graph") which also has minima as the vertices. Gdist has an edge
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
        allow user to pass graph
        
    """    
    def __init__(self, min1, min2, pot, mindist, database, 
                 use_all_min=False, verbosity=1,
                 merge_minima=False, 
                 max_dist_merge=0.1, local_connect_params=dict(),
                 fresh_connect=False, longest_first=False,
                 ):
        self.minstart = min1
        assert min1._id == min1, "minima must compare equal with their id %d %s %s" % (min1._id, str(min1), str(min1.__hash__()))
        self.minend = min2
        self.pot = pot
        self.mindist = mindist
        self.pairsNEB = dict()
        self.longest_first = longest_first
        
        self.verbosity = int(verbosity)
        self.local_connect_params = dict([("verbosity",verbosity)] + local_connect_params.items())
        self.database = database
        self.fresh_connect = fresh_connect
        if self.fresh_connect:
            self.graph = Graph(self.database, minima=[self.minstart, self.minend])
        else:
            self.graph = Graph(self.database)

        self.merge_minima = merge_minima
        self.max_dist_merge = float(max_dist_merge)


        self.dist_graph = _DistanceGraph(self.database, self.graph, self.mindist, self.verbosity)
        self.dist_graph.initialize(self.minstart, self.minend, use_all_min)
        
        print "************************************************************"
        print "starting a double ended connect run between"
        print "        minimum 1: id %d energy %f" % (self.minstart._id, self.minstart.energy)
        print "        minimum 2: id %d energy %f" % (self.minend._id, self.minend.energy)
        print "        dist %f" % self.getDist(self.minstart, self.minend)
        print "************************************************************"
        
    

    
    def mergeMinima(self, min1, min2):
        """merge two minimum objects
        
        This will delete min2 and make everything that
        pointed to min2 point to min1.
        """
        if False:
            print "MERGE MINIMA IS NOT WORKING YET"
            return
        #prefer to delete the minima with the large id.  this potentially will be easier
        if min2._id < min1._id:
            min1, min2 = min2, min1
        
        debug = False
        dist = self.getDist(min1, min2)
        print "merging minima", min1._id, min2._id, dist, "E1-E2", min1.energy - min2.energy

        #deal with the case where min1 and/or min2 are the same as minstart and/or minend
        #make sure the one that is deleted (min2) is not minstart or minend
        if ((min1 == self.minstart and min2 == self.minend) or 
            (min2 == self.minstart and min1 == self.minend)):
            print "ERROR: trying to merge the start and end minima.  aborting"
            return
        if min2 == self.minstart or min2 == self.minend:
            min1, min2 = min2, min1
        
        #print "min1 min2", min1, min2

        if dist > self.max_dist_merge:
            print "    minima merge aborted.  distance is too large", dist
            return
        
#        if debug:
#            #testing
#            if min2 in self.database.minima():
#                print "error, min2 is still in database"
#            for ts in self.database.transition_states():
#                if min2 == ts.minimum1 or min2 == ts.minimum2:
#                    print "error, a transition state attached to min2 is still in database", ts.minimum1._id, ts.minimum2._id
        
        #merge minima in transition state graph
        #note, this will merge minima in the database also
        self.graph.mergeMinima(min1, min2, update_database=True)
        if debug:
            #testing
            if min2 in self.graph.graph.nodes():
                print "error, min2 is still in self.graph.graph"
            print "self.graph.graph. nnodes", self.graph.graph.number_of_nodes()

        #merge minima in distance graph        
        self.dist_graph.mergeMinima(min1, min2)

    def getDist(self, min1, min2):
        """
        get the distance between min1 and min2.
        
        Try first to get distances from the dictionary distmatrix as this is 
        the fastest access method.  Then try to 
        get distances from the database if they exist, else calculate the
        distance and save it to the database and distmatrix
        """
        return self.dist_graph.getDist(min1, min2)

    def _addTransitionState(self, E, coords, min_ret1, min_ret2, eigenvec, eigenval):
        """
        add a transition state to the database, the transition state graph and
        the distance graph
        """
        #sanity check for the energies
        me1, me2 = min_ret1[1], min_ret2[1]
        if E < me1 or E < me2:
            print "warning: trying to add a transition state that has energy lower than it's minima."
            print "    TS energy", E, "minima energy", me1, me2
            print "    aborting"
            return False
        
        #add the minima to the transition state graph.  
        #This step is important to do first because it returns a Database Minimum object.
        min1 = self.graph.addMinimum(min_ret1[1], min_ret1[0])
        min2 = self.graph.addMinimum(min_ret2[1], min_ret2[0])
        if min1 == min2:
            print "warning: stepping off the transition state resulted in twice the same minima", min1._id
            return False
        
        

        print "adding transition state", min1._id, min2._id
        #update the transition state graph
        #this also updates the database and returns a TransitionState object
        ts = self.graph.addTransitionState(E, coords, min1, min2, eigenvec=eigenvec, eigenval=eigenval)
#        self.graph.refresh()

        #update the distance graph
        self.dist_graph.addMinimum(min1)
        self.dist_graph.addMinimum(min2)
        self.dist_graph.setTransitionStateConnection(min1, min2)

        if True:
            #print some information
            dse  = self.getDist(self.minend, self.minstart)
            msid = self.minstart._id
            meid = self.minend._id
            m1id = min1._id
            m2id = min2._id
            if min1 != self.minstart and min1 != self.minend:
                ds = self.getDist(min1, self.minstart)
                de = self.getDist(min1, self.minend)
                if ds < dse > de:
                    triangle = ""
                else: 
                    triangle = ": new minima not in between start and end"
                print "    distances: %4d -> %4d = %f    %4d -> %4d = %f    %4d -> %4d = %f  %s" % (msid, m1id, ds, m1id, meid, de, m1id, m2id, dse, triangle)
                #print "    dist new min 1 to minstart, minend ", ds, de, dse
            if min2 != self.minstart and min2 != self.minend:
                ds = self.getDist(min2, self.minstart)
                de = self.getDist(min2, self.minend)
                if ds < dse > de:
                    triangle = ""
                else: 
                    triangle = ": new minima not in between start and end"
                print "    distances: %4d -> %4d = %f    %4d -> %4d = %f    %4d -> %4d = %f" % (msid, m2id, ds, m2id, meid, de, m2id, m2id, dse)

        
        return True

    def _getLocalConnectObject(self):
        return LocalConnect(self.pot, self.mindist, **self.local_connect_params)

    def _localConnect(self, min1, min2):
        """
        1) NEB to find transition state candidates.  
        
        for each transition state candidate:
        
            2) refine the transition state candidates
        
            3) if successful, fall off either side of the transition state
            to find the minima the transition state connects. Add the new 
            transition state and minima to the graph 
        """
        #Make sure we haven't already tried this pair and
        #record some data so we don't try it again in the future
        if self.pairsNEB.has_key((min1, min2)):
            print "WARNING: redoing NEB for minima", min1._id, min2._id
            print "         aborting NEB"
            #self._remove_edgeGdist(min1, min2)
            self.dist_graph.removeEdge(min1, min2)
            return True
        self.pairsNEB[(min1, min2)] = True
        self.pairsNEB[(min2, min1)] = True
        
        #Make sure they're not already connected.  sanity test
        if self.graph.areConnected(min1, min2):
            print "in _local_connect, but minima are already connected. aborting", min1._id, min2._id, self.getDist(min1, min2)
            self.dist_graph.setTransitionStateConnection(min1, min2)
            self.dist_graph.checkGraph()
            return True
        
        #do local connect run
        local_connect = self._getLocalConnectObject()        
        res = local_connect.connect(min1, min2)

        #now add each new transition state to the graph and database.
        nsuccess = 0
        for tsret, m1ret, m2ret in res.new_transition_states:
            goodts = self._addTransitionState(tsret.energy, tsret.coords, m1ret, m2ret, tsret.eigenvec, tsret.eigenval)
            if goodts:
                nsuccess += 1
 
        #check results
        #nclimbing = len(climbing_images)
        #print "from NEB search found", nclimbing, "transition state candidates"
        if nsuccess == 0:
            dist = self.getDist(min1, min2)
            if dist < self.max_dist_merge:
                print "WARNING: local connect failed and the minima are close. Are the minima really the same?"
                print "         energies:", min1.energy, min2.energy, "distance", dist 
                if self.merge_minima:
                    self.mergeMinima(min1, min2)
                else:
                    print "         set merge_minima=True to merge the minima" 
                return False   


        #remove this edge from Gdist so we don't try this pair again again
        #self._remove_edgeGdist(min1, min2)
        self.dist_graph.removeEdge(min1, min2)
        
        return nsuccess > 0            
    
                
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
        #get the shortest path on dist_graph between minstart and minend
        path, weights = self.dist_graph.shortestPath(self.minstart, self.minend)
        weightsum = sum(weights)
        if path is None or weightsum >= 10e9:
            print "Can't find any way to try to connect the minima"
            return None, None
        
        #get the weights of the path segements
        weightlist = []
        for i in range(1,len(path)):
            min1 = path[i-1]
            min2 = path[i]
#            w = weights.get((min1,min2))
#            if w is None:
#                w = weights.get((min2,min1))
            w = weights[i-1]
            weightlist.append( (w, min1, min2) )
        
        if True:
            #print the path
            print "best guess for path.  (dist=0.0 means the path is known)"
            for w, min1, min2 in weightlist:
                if w > 1e-6:
                    dist = self.getDist(min1, min2)
                else:
                    dist = w
                print "    path guess", min1._id, min2._id, dist

        #select which minima pair to return
        if self.longest_first:
            weightlist.sort()
            w, min1, min2 = weightlist[-1]
        else:
            weightlist.sort()
            for w, min1, min2 in weightlist:
                if w > 1e-6:
                    break
        return min1, min2
            
        
#        #find the shortest single edge in that path and use that
#        #also, print the path
#        #as the minima pair to try to connect.
#        weightmin = 1e100
#        minpair = (None, None)
#        for i in range(1,len(path)):
#            min1 = path[i-1]
#            min2 = path[i]
#            w = weights.get((min1,min2))
#            if w is None:
#                w = weights.get((min2,min1))
#            #get the distance between min1 and min2. leave as 0. if they are
#            #already connected by transition states
#            dist = w
#            if w > 1e-6:
#                dist = self.getDist(min1, min2)
#            print "    path guess", min1._id, min2._id, dist
#            if w > 1e-10 and w < weightmin:
#                weightmin = w
#                minpair = (min1, min2)
#        return minpair
                
    
    def connect(self, maxiter=200):
        """
        the main loop of the algorithm
        """
        self.NEBattempts = 2;
        for i in range(maxiter):
            self.dist_graph.updateDatabase()
            if self.graph.areConnected(self.minstart, self.minend):
                self.dist_graph.updateDatabase(force=True)
                print "found connection!"
                return
            
            print ""
            print "======== starting connect cycle", i, "========"
            min1, min2 = self._getNextPair()
            #raw_input("Press Enter to continue:")
            if min1 is None or min2 is None:
                break
            local_success = self._localConnect(min1, min2)
            
            if True and i % 10 == 0:
                #do some santy checks
                self.dist_graph.checkGraph()



            
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

def getSetOfMinLJ(natoms = 32): #for testing purposes
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


def test(Connect=DoubleEndedConnect, natoms=16):
    from pygmin.landscape import Graph
    from pygmin.optimize.quench import lbfgs_py as quench
    from pygmin.mindist.minpermdist_stochastic import minPermDistStochastic as mindist
    from pygmin.storage.database import Database
    import pygmin.defaults as defaults
    defaults.quenchParams = {"iprint": 1}
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
 
    connect = Connect(min1, min2, pot, mindist, database)
    connect.connect()
    
    graph = connect.graph
    if False:
        print graph
        for node in graph.graph.nodes():
            print node._id, node.energy
    for ts in graph.storage.transition_states():
        print ts.minimum1._id,ts.minimum2._id, "E", ts.minimum1.energy, ts.minimum2.energy, ts.energy
        
    ret = graph.getPath(min1, min2)
    if ret is None:
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
    test(natoms=38)


    
