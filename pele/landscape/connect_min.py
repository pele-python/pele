from __future__ import print_function
import logging
import operator

import numpy as np
import networkx as nx

from pele.landscape import TSGraph, LocalConnect
from pele.landscape._distance_graph import _DistanceGraph


__all__ = ["DoubleEndedConnect"]

logger = logging.getLogger("pele.connect")


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
    database : pele Database object
        Used to store the new minima and transition states found.
    niter : int, optional
        maximum number of iterations
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
    conf_checks : list of callables
        a list of callable function that determine if a configuration is valid.
        They must return a bool, and accept the keyword parameters
        
            conf_check(energy=energy, coords=coords)
        
        If any configuration in a minimum-transition_state-minimum triplet fails
        a test then the whole triplet is rejected.
    
    Notes
    -----
    The algorithm is iterative, with each iteration composed of
    
    While min1 and min2 are not connected:
        1) choose a pair of known minima to try to connect
        
        2) use NEB to get a guess for the transition states between them
        
        3) refine the transition states to desired accuracy
        
        4) fall off either side of the transition states to find the two
        minima associated with that candidate
        
        5) add the transition states and associated minima to the known
        network
        
    
    Of the above, steps 1 and 2 and 3 are the most involved.  2, 3, 4 are 
    wrapped into a separate class called LocalConnect.  See this class and
    the NEB and FindTransitionState classes for detailed descriptions of 
    these procedures.
    
    An important note is that the NEB is used only to get a *guess* for the
    transition state.  Thus we only want to put enough time and energy into
    the NEB routine to get the guess close enough that FindTransitionState
    can refine it to the correct transition state.  FindTransitionState is
    very fast if the initial guess is good, but can be very slow otherwise.
    
    Choose a pair:

    Here I will describe step 1), the algorithm to find a pair of known
    minima to try to connect.  This choice will keep in mind that the
    ultimate goal is to connect min1 and min2.
    
    In addition to the input parameter "graph", we keep a second graph
    "Gdist" (now wrapped in a separate class _DistanceGraph) which also has 
    minima as the vertices. Gdist has an edge between every pair of nodes. 
    The edge weight between vertices u and v
    is
    
        if u and v are connected by transition states:
            weight(u, v) = 0. 
        elif we have already tried local_connect on (u,v):
            weight(u, v) = Infinity
        else:
            weight(u, v) = dist(u, v)**2
    
    This edge weight is set to Infinity to ensure we don't repeat 
    LocalConnect runs over and over
    again.  The minimum weight path between min1 and min2 in Gdist gives a
    good guess for the best way to try connect min1 and min2.  So the
    algorithm to find a pair of know minima (trial1, trial2) to try to
    connect is 
    
    path = Gdist.minimum_weight_path(min1, min2)
    trial1, trial2 = minima pair in path with lowest nonzero edge weight. (note:
    if parameter longest_first is True) then the edgepair with the largest
    edge weight will be selected) 

    

    todo:
        allow user to pass graph
    
    See Also
    --------
    LocalConnect : the core algorithm of this routine
        
    """

    def __init__(self, min1, min2, pot, mindist, database,
                 verbosity=1,
                 merge_minima=False,
                 max_dist_merge=0.1, local_connect_params=None,
                 fresh_connect=False, longest_first=True,
                 niter=200, conf_checks=None
    ):
        self.minstart = min1
        assert min1.id() == min1, "minima must compare equal with their id %d %s %s" % (
            min1.id(), str(min1), str(min1.__hash__()))
        self.minend = min2
        self.pot = pot
        self.mindist = mindist
        self.pairsNEB = dict()
        self.longest_first = longest_first
        self.niter = niter
        if conf_checks is None:
            self.conf_checks = []
        else:
            self.conf_checks = conf_checks

        self.verbosity = int(verbosity)
        if local_connect_params is None:
            local_connect_params = dict()
        self.local_connect_params = dict([("verbosity", verbosity)] + list(local_connect_params.items()))
        self.database = database
        self.fresh_connect = fresh_connect
        if self.fresh_connect:
            self.graph = TSGraph(self.database, minima=[self.minstart, self.minend], no_edges=True)
        else:
            self.graph = TSGraph(self.database)

        self.merge_minima = merge_minima
        self.max_dist_merge = float(max_dist_merge)

        self.dist_graph = _DistanceGraph(self.database, self.graph, self.mindist, self.verbosity)

        # check if a connection exists before initializing distance graph
        if self.graph.areConnected(self.minstart, self.minend):
            logger.info("minima are already connected.  not initializing distance graph")
            return

        self.dist_graph.initialize(self.minstart, self.minend)

        if self.verbosity > 0:
            logger.info("************************************************************")
            logger.info("starting a double ended connect run between")
            logger.info("        minimum 1: id %d energy %f" % (self.minstart.id(), self.minstart.energy))
            logger.info("        minimum 2: id %d energy %f" % (self.minend.id(), self.minend.energy))
            logger.info("        dist %f" % self.getDist(self.minstart, self.minend))
            logger.info("************************************************************")


    def mergeMinima(self, min1, min2):
        """merge two minimum objects
        
        This will delete min2 and make everything that
        pointed to min2 point to min1.
        """
        # prefer to delete the minima with the large id.  this potentially will be easier
        if min2.id() < min1.id():
            min1, min2 = min2, min1

        debug = False
        dist = self.getDist(min1, min2)
        logger.info("merging minima %s %s %s %s %s", min1.id(), min2.id(), dist, "E1-E2", min1.energy - min2.energy)

        # deal with the case where min1 and/or min2 are the same as minstart and/or minend
        # make sure the one that is deleted (min2) is not minstart or minend
        if ((min1 == self.minstart and min2 == self.minend) or
                (min2 == self.minstart and min1 == self.minend)):
            logger.error("ERROR: trying to merge the start and end minima.  aborting")
            return
        if min2 == self.minstart or min2 == self.minend:
            min1, min2 = min2, min1

        if dist > self.max_dist_merge:
            logger.info("    minima merge aborted.  distance is too large %s", dist)
            return

        # merge minima in transition state graph
        self.graph.mergeMinima(min1, min2)

        # merge minima in the database also
        self.database.mergeMinima(min1, min2)
        if debug:
            # testing
            if min2 in self.graph.graph.nodes():
                logger.error("error, min2 is still in self.graph.graph")
            logger.debug("self.graph.graph. nnodes %s", self.graph.graph.number_of_nodes())

        # merge minima in distance graph        
        self.dist_graph.mergeMinima(min1, min2)

    def getDist(self, min1, min2):
        """
        get the distance between min1 and min2.
        """
        return self.dist_graph.getDist(min1, min2)

    def _addTransitionState(self, ts_ret, min_ret1, min_ret2):
        """
        add a transition state to the database, the transition state graph and
        the distance graph
        """
        # if isinstance(min_ret1, tuple): # for compatability with old and new quenchers
        # min_ret1 = min_ret1[4]
        # if isinstance(min_ret2, tuple): # for compatability with old and new quenchers
        # min_ret2 = min_ret2[4]
        # if isinstance(min1_ret)

        # sanity check for the energies
        me1, me2 = min_ret1.energy, min_ret2.energy
        if ts_ret.energy < me1 or ts_ret.energy < me2:
            logger.warning("trying to add a transition state that has energy lower than its minima.")
            logger.warning("    TS energy %s %s %s %s", ts_ret.energy, "minima energy", me1, me2)
            logger.warning("    aborting")
            return False

        # check the minima and transition states are valid configurations.
        # if any fail, then don't add anything.  
        configs_ok = True
        for ret in [min_ret1, min_ret2, ts_ret]:
            for check in self.conf_checks:
                if not check(energy=ret.energy, coords=ret.coords):
                    configs_ok = False
                    break
            if not configs_ok:
                break
        if not configs_ok:
            return False


        # Add the minima to the database
        min1 = self.database.addMinimum(min_ret1.energy, min_ret1.coords)
        min2 = self.database.addMinimum(min_ret2.energy, min_ret2.coords)

        # Add the minima to the transition state graph.  
        self.graph.addMinimum(min1)
        self.graph.addMinimum(min2)
        if min1 == min2:
            logger.warning("stepping off the transition state resulted in twice the same minima %s", min1.id())
            return False

        logger.info("adding transition state %s %s", min1.id(), min2.id())
        # add the transition state to the database
        ts = self.database.addTransitionState(ts_ret.energy, ts_ret.coords, min1, min2,
                                              eigenval=ts_ret.eigenval)
        # update the transition state graph
        self.graph.addTransitionState(ts)
        # self.graph.refresh()

        # update the distance graph
        self.dist_graph.addMinimum(min1)
        self.dist_graph.addMinimum(min2)
        self.dist_graph.setTransitionStateConnection(min1, min2)

        if self.verbosity > 1:
            # print some information
            dse = self.getDist(self.minend, self.minstart)
            msid = self.minstart.id()
            meid = self.minend.id()
            m1id = min1.id()
            m2id = min2.id()
            if min1 != self.minstart and min1 != self.minend:
                ds = self.getDist(min1, self.minstart)
                de = self.getDist(min1, self.minend)
                if ds < dse > de:
                    triangle = ""
                else:
                    triangle = ": new minima not in between start and end"
                logger.info("    distances: %4d -> %4d = %f    %4d -> %4d = %f    %4d -> %4d = %f  %s" %
                            (msid, m1id, ds, m1id, meid, de, m1id, m2id, dse, triangle))
            if min2 != self.minstart and min2 != self.minend:
                ds = self.getDist(min2, self.minstart)
                de = self.getDist(min2, self.minend)
                # if ds < dse > de:
                # triangle = ""
                # else:
                # triangle = ": new minima not in between start and end"
                logger.info("    distances: %4d -> %4d = %f    %4d -> %4d = %f    %4d -> %4d = %f" %
                            (msid, m2id, ds, m2id, meid, de, m2id, m2id, dse))

        return True

    def _getLocalConnectObject(self):
        return LocalConnect(self.pot, self.mindist, **self.local_connect_params)

    def _localConnect(self, min1, min2):
        """
        do a local connect run between min1 and min2
        
        Notes
        -----
        1) NEB to find transition state candidates.  
        
        for each transition state candidate:
        
            2) refine the transition state candidates
        
            3) if successful, fall off either side of the transition state
            to find the minima the transition state connects. Add the new 
            transition state and minima to the graph 
        """
        # Make sure we haven't already tried this pair and
        # record some data so we don't try it again in the future
        if (min1, min2) in self.pairsNEB:
            logger.warning("WARNING: redoing NEB for minima %s %s", min1.id(), min2.id())
            logger.warning("         aborting NEB")
            # self._remove_edgeGdist(min1, min2)
            self.dist_graph.removeEdge(min1, min2)
            return True
        self.pairsNEB[(min1, min2)] = True
        self.pairsNEB[(min2, min1)] = True

        # Make sure they're not already connected.  sanity test
        if self.graph.areConnected(min1, min2):
            logger.warning("in _local_connect, but minima are already connected. aborting %s %s %s", min1.id(), min2.id(),
                           self.getDist(min1, min2))
            self.dist_graph.setTransitionStateConnection(min1, min2)
            self.dist_graph.checkGraph()
            return True

        # do local connect run
        local_connect = self._getLocalConnectObject()
        res = local_connect.connect(min1, min2)

        # now add each new transition state to the graph and database.
        nsuccess = 0
        for tsret, m1ret, m2ret in res.new_transition_states:
            goodts = self._addTransitionState(tsret, m1ret, m2ret)
            if goodts:
                nsuccess += 1

        # check results
        if nsuccess == 0:
            dist = self.getDist(min1, min2)
            if dist < self.max_dist_merge:
                logger.warning("local connect failed and the minima are close. Are the minima really the same?")
                logger.warning("         energies: %s %s %s %s", min1.energy, min2.energy, "distance", dist)
                if self.merge_minima:
                    self.mergeMinima(min1, min2)
                else:
                    logger.warning("         set merge_minima=True to merge the minima")
                return False


                # remove this edge from Gdist so we don't try this pair again again
        self.dist_graph.removeEdge(min1, min2)

        return nsuccess > 0


    def _getNextPair(self):
        """
        return a pair of minima to attempt to connect
        
        Notes
        -----
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
        logger.info("finding a good pair to try to connect")
        # get the shortest path on dist_graph between minstart and minend
        if True:
            logger.debug("Gdist has %s %s %s %s", self.dist_graph.Gdist.number_of_nodes(),
                         "nodes and", self.dist_graph.Gdist.number_of_edges(), "edges")
        path, weights = self.dist_graph.shortestPath(self.minstart, self.minend)
        weightsum = sum(weights)
        if path is None or weightsum >= 10e9:
            logger.warning("Can't find any way to try to connect the minima")
            return None, None

        # get the weights of the path segements
        weightlist = []
        for i in range(1, len(path)):
            min1 = path[i - 1]
            min2 = path[i]
            # w = weights.get((min1,min2))
            # if w is None:
            # w = weights.get((min2,min1))
            w = weights[i - 1]
            weightlist.append((w, min1, min2))

        if True:
            # print the path
            logger.info("best guess for path.  (dist=0.0 means the path is known)")
            for w, min1, min2 in weightlist:
                if w > 1e-6:
                    dist = self.getDist(min1, min2)
                else:
                    dist = w
                logger.info("    path guess %s %s %s", min1.id(), min2.id(), dist)

        # select which minima pair to return
        if self.longest_first:
            weightlist.sort(key=operator.itemgetter(0))
            w, min1, min2 = weightlist[-1]
        else:
            weightlist.sort(key=operator.itemgetter(0))
            for w, min1, min2 in weightlist:
                if w > 1e-6:
                    break
        return min1, min2


    def connect(self):
        """
        the main loop of the algorithm
        """
        self.NEBattempts = 2
        for i in range(self.niter):
            # stop if we're done
            if self.graph.areConnected(self.minstart, self.minend):
                logger.info("found connection!")
                return

            logger.info("")
            logger.info("======== starting connect cycle %s %s", i, "========")
            # get pair of minima to try to connect
            min1, min2 = self._getNextPair()

            # fail if we can't find a good pair to try
            if min1 is None or min2 is None:
                break

            # try to connect those minima
            from pele.optimize.optimization_exceptions import LineSearchError

            try:
                self._localConnect(min1, min2)
            except LineSearchError as err:
                print(err)
                print("caught line search error, aborting connection attempt")
                break

            if False and i % 10 == 0:
                # do some sanity checks
                self.dist_graph.checkGraph()

        logger.info("failed to find connection between %s %s", self.minstart.id(), self.minend.id())

    def success(self):
        return self.graph.areConnected(self.minstart, self.minend)

    def returnPath(self):
        """return information about the path
        
        Returns
        -------
        mints : list of Minimum and TransitionStates
            a list of Minimum, TransitionState, Minimum objects that make up
            the path
        S : list of float 
            numpy array of the distance along the path.   len(S) == len(mints)
        energies : list of float
            numpy array of the energies along the path

        If the minima are not connected, return (None, None, None)
        """
        if not self.graph.areConnected(self.minstart, self.minend):
            return None, None, None
        minima = nx.shortest_path(self.graph.graph, self.minstart, self.minend)
        transition_states = []
        mints = [minima[0]]
        for i in range(1, len(minima)):
            m1 = minima[i - 1]
            m2 = minima[i]
            ts = self.database.getTransitionState(m1, m2)
            transition_states.append(ts)
            mints.append(ts)
            mints.append(m2)

        S = np.zeros(len(mints))
        for i in range(1, len(mints)):
            coords1 = mints[i - 1].coords
            coords2 = mints[i].coords
            dist, c1, c2 = self.mindist(coords1, coords2)
            S[i] = S[i - 1] + dist
        energies = np.array([m.energy for m in mints])
        return mints, S, energies


# ##########################################################
# only testing stuff below here
# ##########################################################

def getSetOfMinLJ(system):  # pragma: no cover
    db = system.create_database()
    bh = system.get_basinhopping(db, outstream=None)
    bh.run(100)
    return system.get_potential(), db


def test(Connect=DoubleEndedConnect, natoms=16):  # pragma: no cover
    from pele.systems import LJCluster
    # get min1
    system = LJCluster(natoms)
    pot, database = getSetOfMinLJ(system)
    minima = database.minima()
    min1 = minima[0]
    min2 = minima[1]
    print(min1.energy, min2.energy)

    mindist = system.get_mindist()

    connect = Connect(min1, min2, pot, mindist, database)
    connect.connect()

    graph = connect.graph
    if False:
        print(graph)
        for node in graph.graph.nodes():
            print(node.id(), node.energy)
    for ts in graph.storage.transition_states():
        print(ts.minimum1.id(), ts.minimum2.id(), "E", ts.minimum1.energy, ts.minimum2.energy, ts.energy)

    ret = graph.getPath(min1, min2)
    if ret is None:
        print("no path found")
        return
    distances, path = ret
    with open("path.out", "w") as fout:
        for i in range(len(path) - 1):
            m1 = path[i]
            m2 = path[i + 1]
            n1 = m1.id()
            m2 = m2.id()
            # ts = graph._getTS(n1, n2)
            # print "path", n1, "->", n2, m1.E, "/->", ts.E, "\->", m2.E
            fout.write("%f\n" % m1.energy)
            fout.write("%f\n" % ts.energy)
        m2 = path[-1]
        n2 = m2.id()
        fout.write("%f\n" % m2.energy)


if __name__ == "__main__":
    # logger.basicConfig(level=logger.DEBUG)
    test(natoms=38)

    

