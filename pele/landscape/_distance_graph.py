import networkx as nx
import logging


__all__ = []

logger = logging.getLogger("pele.connect")


class _DistanceGraph(object):
    """
    This graph is used by DoubleEndedConnect to make educated guesses for connecting two minima
    
    Parameters
    ----------
    database : 
        the database in which to store newly found minima and transition states.
    graph :
        the graph build from the database which contains the minima and transition states
    mindist : callable
        the routine which calculates the optimized distance between two structures
    verbosity :
        how much info to print (not very thoroughly implemented)
    
    Description
    -----------
    This graph has a vertex for every minimum and an edge
    between every minima pair. The edge weight between vertices u and v
    is
                       
        if u and v are connected by transition states:
            weight(u, v) = 0. 
        elif we have already tried local_connect on (u,v):
            weight(u, v) = Infinity
        else:
            weight(u, v) = dist(u, v)**2

    Edge weights are set to Infinity to ensure that we don't try to connect
    them again.  The minimum weight path between min1 and min2 in this graph gives a
    good guess for the best way to try connect min1 and min2.  

    TODO: could we use a disjoint set data structure (nx.utils.UnionFind) to improve this
    algorithm?
    """

    def __init__(self, database, graph, mindist, verbosity):
        self.database = database
        self.graph = graph
        self.mindist = mindist
        self.verbosity = verbosity

        self.Gdist = nx.Graph()
        self.distance_map = dict()  # place to store distances locally for faster lookup
        nx.set_edge_attributes(self.Gdist, "weight", dict())
        self.debug = False

        self.infinite_weight = 1e20

    def distToWeight(self, dist):
        """
        this function defines how the edge weight is calculated.
        
        good options might be:
        
        weight = dist
        
        weight = dist**2
            this favors paths with many short edges over
            paths with fewer but longer edges.
        """
        return dist ** 2

    def _setDist(self, min1, min2, dist):
        """
        this function saves newly calculated distances both to the local
        distance map
        """
        # add the distance to the distance map
        self.distance_map[(min1, min2)] = dist

        # make sure a zeroed edge weight is not overwritten
        # if not self.edge_weight.has_key((min1, min2)):
        # if not self.edge_weight.has_key((min2, min1)):
        # self.edge_weight[(min1, min2)] = weight

    def _getDistNoCalc(self, min1, min2):
        """
        get distance from local memory.  if it doesn't exist, return None,
        don't calculate it.
        """
        # first try to get the distance from the dictionary 
        dist = self.distance_map.get((min1, min2))
        if dist is not None: return dist
        dist = self.distance_map.get((min2, min1))
        if dist is not None: return dist
        return None

    def getDist(self, min1, min2):
        """
        return the distance between two minima.  Calculate it and store it if
        not already known
        """
        dist = self._getDistNoCalc(min1, min2)
        if dist is not None: return dist

        # if it's not already known we must calculate it
        dist, coords1, coords2 = self.mindist(min1.coords, min2.coords)
        if self.verbosity > 1:
            logger.debug("calculated distance between %s %s %s", min1.id(), min2.id(), dist)
        self._setDist(min1, min2, dist)
        return dist

    def _addMinimum(self, m):
        """
        add a new minimum to the graph
        
        must add an edge with the appropriate weight to every other 
        node in the graph.
        
        this can take a long time if there are many minima or if the
        distance calculation is slow.
        """
        self.Gdist.add_node(m)
        # for nodes that are connected set the edge weight using setTransitionStateConnection
        cc = nx.node_connected_component(self.graph.graph, m)
        for m2 in cc:
            if m2 in self.Gdist:
                # self.Gdist.add_edge(m, m2, weight=0.)
                self.setTransitionStateConnection(m, m2)

        # for all other nodes set the weight to be the distance
        for m2 in self.Gdist.nodes():
            if not self.Gdist.has_edge(m, m2):
                dist = self.getDist(m, m2)
                weight = self.distToWeight(dist)
                self.Gdist.add_edge(m, m2, weight = weight)


    def addMinimum(self, m):
        """
        add a new minima to the graph and add an edge to all the other
        minima in the graph.  
        
        Note: this can take a very long time if there are lots of minima
        in the graph.  mindist need to be run many many times.
        """
        trans = self.database.connection.begin()
        try:
            if m not in self.Gdist:
                self._addMinimum(m)
        except BaseException:
            trans.rollback()
            raise
        trans.commit()


    def removeEdge(self, min1, min2):
        """set the edge weight to near infinity
                
        used to indicate that the routine should not try to connect
        these minima again.
        
        don't overwrite zero edge weight
        """
        if self.Gdist.has_edge(min1, min2):
            w = self.Gdist[min1][min2]["weight"]
            if not w < 1e-6:
                self.Gdist.add_edge(min1, min2, weight=self.infinite_weight)
        return True

    def replaceTransitionStateGraph(self, graph):
        self.graph = graph

    def initialize(self, minstart, minend):
        """
        set up the distance graph
        
        initialize distance_map, add the start and end minima and load any other
        minima that should be used in the connect routine.
        """
        self.getDist(minstart, minend)
        self.addMinimum(minstart)
        self.addMinimum(minend)

    def setTransitionStateConnection(self, min1, min2):
        """use this function to tell _DistanceGraph that
        there exists a known transition state connection between min1 and min2
        
        The edge weight will be set to zero
        """
        weight = 0.
        self.Gdist.add_edge(min1, min2, weight = weight)

    def shortestPath(self, min1, min2):
        """return the minimum weight path path between min1 and min2"""
        try:
            path = nx.shortest_path(self.Gdist, min1, min2, weight="weight")
        except nx.NetworkXNoPath:
            return None, None

        # get_edge attributes is really slow:
        weights = [self.Gdist[path[i]][path[i + 1]]["weight"] for i in range(len(path) - 1)]

        return path, weights

    def mergeMinima(self, min1, min2):
        """
        rebuild the graph with min2 deleted and 
        everything pointing to min1 pointing to min2 instead
        """
        for m, data in self.Gdist[min2].items():
            if m == min1:
                continue
            # if not self.Gdist.has_edge(min1, m):
            # self.add_edge(min1, m, **data)

            # the edge already exists, keep the edge with the lower weight
            w2 = data["weight"]
            w1 = self.Gdist[min1][m]["weight"]
            wnew = min(w1, w2)
            # note: this will override any previous call to self.setTransitionStateConnection
            self.Gdist.add_edge(min1, m, weight=wnew)

        self.Gdist.remove_node(min2)


    def checkGraph(self):
        """
        make sure graph is up to date and make any corrections
        """
        logger.info("checking Gdist")
        allok = True
        # check that all edges that are connected in self.graph
        # have zero edge weight
        # note: this could be done a lot more efficiently
        weights = nx.get_edge_attributes(self.Gdist, "weight")
        count = 0
        for e in self.Gdist.edges():
            are_connected = self.graph.areConnected(e[0], e[1])
            zero_weight = weights[e] < 1e-10

            # if they are connected they should have zero_weight
            if are_connected and not zero_weight:
                # print "    problem: are_connected", are_connected, "but weight", weights[e], "dist", dist
                if True:
                    # note: this is an inconsistency, but it's only a problem if
                    # there is no zero weight path from e[0] to e[1]
                    path, path_weight = self.shortestPath(e[0], e[1])
                    weight_sum = sum(path_weight)
                    if weight_sum > 10e-6:
                        # now there is definitely a problem.
                        allok = False
                        count += 1
                        dist = self.getDist(e[0], e[1])
                        logger.warning("    problem: are_connected %s %s %s %s %s %s %s %s %s",
                                       are_connected, "but weight", weights[e], "dist", dist,
                                       "path_weight", weight_sum, e[0].id(), e[1].id())
                self.setTransitionStateConnection(e[0], e[1])

            if not are_connected and zero_weight:
                allok = False
                dist = self.getDist(e[0], e[1])
                logger.warning("    problem: are_connected %s %s %s %s %s %s %s",
                               are_connected, "but weight", weights[e], "dist", dist, e[0].id(), e[1].id())
                w = self.distToWeight(dist)
                self.Gdist.add_edge(e[0], e[1], weight = w)
        if count > 0:
            logger.info("    found %s %s", count, "inconsistencies in Gdist")

        return allok

