"""
classes to organize strategies for selecting which minima in a database to
choose for a double ended connect run. 
"""
from __future__ import print_function
from collections import deque

import numpy as np
import sqlalchemy
import networkx as nx

from pele.storage import Minimum
from pele.landscape import TSGraph


__all__ = ["ConnectManager"]


class BaseConnectManager(object):
    _is_good_pair = lambda self, m1, m2: True

    list_len = 10
    minpairs = None

    """
    Attributes
    ----------
    list_len : int
        the class will create a list of minima pairs of length
        list_len.  When this list is empty the list will be rebuilt.
        Essentially this parameter indicates how often to rebuild the list
        of minima pairs to connect
    minpairs : deque
        This will is the queue of minima pairs to try
    """

    def _build_list(self):
        return NotImplementedError

    def set_good_pair_test(self, test):
        self._is_good_pair = test

    def is_good_pair(self, min1, min2):
        return self._is_good_pair(min1, min2)

    def get_connect_job(self):
        if len(self.minpairs) == 0:
            self._build_list()
        if len(self.minpairs) == 0:
            # this method has run out of options.  The next list_len jobs will be (None, None)
            self.minpairs = deque([(None, None)] * self.list_len)
            return None, None

        min1, min2 = self.minpairs.popleft()
        return min1, min2


class ConnectManagerGMin(BaseConnectManager):
    """
    return connect jobs in order to connect everything to the global minimum
    """

    def __init__(self, database, list_len=10, verbosity=1):
        self.database = database
        self.list_len = list_len
        self.verbosity = verbosity

        self.minpairs = deque()

    def _build_list(self):
        if self.verbosity > 0:
            print("populating list of minima not connected to the global minimum")
        self.minpairs = deque()

        gmin = self.database.minima()[0]

        graph = TSGraph(self.database)

        for m in self.database.minima()[1:]:
            if not graph.areConnected(gmin, m):
                if self.is_good_pair(gmin, m):
                    self.minpairs.append((gmin, m))


class ConnectManagerUntrap(BaseConnectManager):
    """class to select double ended connect jobs using the untrap strategy
    
    Notes
    ------
    the untrap strategy selects minima with the goal of eliminating artificially high
    energy barriers.  Each minima is weighted by the energy barrier to get to the global
    minimum divided by the energy difference with the global minimum.  This weight determines
    the order with in the minima are selected for a double ended connect run.
    
    Parameters
    ----------
    database : pele Database object
    list_len : int
        the class will create a list of minima pairs of length
        list_len.  When this list is empty the list will be rebuilt.
        Essentially this parameter indicates how often to rebuild the list
        of minima pairs to connect
    """

    def __init__(self, database, list_len=10, nlevels=20):
        self.database = database
        self.list_len = list_len
        self.nlevels = nlevels

        self.minpairs = deque()

    def _recursive_label(self, tree, min1, energy_barriers):
        if tree.is_leaf():
            return

        for subtree in tree.subtrees:
            if subtree.contains_minimum(min1):
                self._recursive_label(subtree, min1, energy_barriers)
            else:
                energy_barrier = tree.data["ethresh"]
                for min2 in subtree.get_minima():
                    assert min2 != min1
                    # print "minimum", min2.id(), "has energy barrier", energy_barrier
                    energy_barriers[min2] = energy_barrier - min2.energy

    def _compute_barriers(self, graph, min1):
        """for each minimum graph compute the (approximate) energy barrier to min1"""
        # this is a local import to avoid cyclical imports
        from pele.utils.disconnectivity_graph import DisconnectivityGraph

        dgraph = DisconnectivityGraph(graph, nlevels=self.nlevels)
        dgraph.calculate()
        tree = dgraph.tree_graph

        energy_barriers = dict()
        self._recursive_label(tree, min1, energy_barriers)
        return energy_barriers

    def _build_list(self):
        print("using disconnectivity analysis to find minima to untrap")
        self.minpairs = deque()

        graph = TSGraph(self.database).graph
        cclist = list(nx.connected_components(graph))

        # get the largest cluster
        group1 = cclist[0]
        min1 = sorted(group1, key=lambda m: m.energy)[0]
        if not min1 == self.database.minima()[0]:
            # make sure that the global minimum is in group1
            print("warning, the global minimum is not the in the largest cluster.")

        # compute the energy barriers for all minima in the cluster        
        subgraph = nx.subgraph(graph, group1)
        energy_barriers = self._compute_barriers(subgraph, min1)

        # sort the minima by the barrier height divided by the energy difference
        weights = [(m, np.abs(barrier) / np.abs(m.energy - min1.energy))
                   for (m, barrier) in energy_barriers.items()]
        weights.sort(key=lambda v: 1. / v[1])

        self.minpairs = deque()
        for min2, w in weights:
            if len(self.minpairs) > self.list_len:
                break

            if not self.is_good_pair(min1, min2):
                continue

            self.minpairs.append((min1, min2))
            if True:
                # print some stuff
                print("    untrap analysis: minimum", min2.id(), "with energy", min2.energy, "barrier", energy_barriers[
                    min2], "untrap weight", w)


class ConnectManagerCombine(BaseConnectManager):
    """a class to organize choosing minima in order to combine disconnected clusters of minima
    
    Parameters
    ----------
    database : Database object
    list_len : int
        the class will create a list of minima pairs of length
        list_len.  When this list is empty the list will be rebuilt.
        Essentially this parameter indicates how often to rebuild the list
        of minima pairs to connect
    clust_min : int
        Clusters of minima below this size will be ignored.
    """

    def __init__(self, database, list_len=20, clust_min=4):
        self.database = database
        self.list_len = list_len
        self.clust_min = clust_min

        self.minpairs = deque()

    def _build_list(self):
        """make a list of minima pairs to try to connect"""
        print("analyzing the database to find minima to connect")
        self.minpairs = deque()

        graph = TSGraph(self.database).graph
        cclist = list(nx.connected_components(graph))

        # remove clusters with fewer than clust_min
        cclist = [cc for cc in cclist if len(cc) >= self.clust_min]

        if len(cclist) == 0:
            print("all minima are connected")
            return self.minpairs

        # get the group that all other groups will be connected to
        group1 = cclist[0]
        min1 = sorted(group1, key=lambda m: m.energy)[0]
        if True:
            # make sure that the global minimum is in group1
            global_min = self.database.minima()[0]
            if global_min not in group1:
                print("warning, the global minimum is not the in the largest cluster.  Will try to connect them")
                self.minpairs.append((min1, global_min))


        # remove group1 from cclist
        cclist.remove(group1)

        # get a minimum from each of the other groups
        for group2 in cclist:
            if len(self.minpairs) > self.list_len:
                break

            print("adding groups of size", len(group1), "and", len(group2), "to the connect list")

            # sort the groups by energy
            group2 = sorted(group2, key=lambda m: m.energy)

            # select the lowest energy minima in the groups
            # (this can probably be done in a more intelligent way)
            min2 = group2[0]

            if self.is_good_pair(min1, min2):
                self.minpairs.append((min1, min2))

        return self.minpairs


class ConnectManagerRandom(BaseConnectManager):
    """manager to return random minima to connect"""

    def __init__(self, database, Emax=None):
        self.database = database
        self.Emax = Emax

    def get_connect_job(self):
        """select two minima randomly"""
        query = self.database.session.query(Minimum)
        if self.Emax is not None:
            query.filter(Minimum.energy < self.Emax)

        itermax = self.database.number_of_minima()  # so we don't iterate forever
        for i in range(itermax):
            min1 = query.order_by(sqlalchemy.func.random()).first()
            min2 = query.order_by(sqlalchemy.func.random()).first()

            if min1 == min2: continue
            if self.is_good_pair(min1, min2):
                return min1, min2

            #        print "worker requested new job, sending minima", min1.id(), min2.id()

        print("warning: couldn't find any random minima pair to connect")
        return None, None


class ConnectManager(object):
    """class to manage which minima to try to connect
    
    Notes
    -----
    Organize which minima pairs to submit for double ended connect jobs.
    This class chooses between the different selection strategies.  
    The actual strategies are implemented in separate classes.

    Selection strategies:
    
        1. "random" : choose two minima randomly
        #. "gmin" : connect all minima to the global minimum
        #. "combine" : try to connect disconnected clusters of minima.
        #. "untrap" : try to eliminate non-physical barriers to the global minimum 
    
    
    Parameters
    ----------
    database : pele Database object
    strategy : string
        define the default strategy for the connect runs.  Can be one of 
        ["random", "combine", "untrap", "gmin"] 
    """

    class NoMoreConnectionsError(Exception):
        """raised when the connect manager can't find any more pairs to connect"""

    def __init__(self, database, strategy="random", list_len=20, clust_min=4, Emax=None,
                 untrap_nlevels=20, verbosity=1):
        self.database = database
        self.default_strategy = strategy
        self.verbosity = verbosity

        self.manager_random = ConnectManagerRandom(self.database, Emax)
        self.manager_combine = ConnectManagerCombine(self.database, list_len=list_len, clust_min=4)
        self.manager_untrap = ConnectManagerUntrap(database, list_len=list_len, nlevels=untrap_nlevels)
        self.manager_gmin = ConnectManagerGMin(database, list_len=list_len, verbosity=self.verbosity)

        self.possible_strategies = ["random", "combine", "untrap", "gmin"]
        self.backup_strategy = "random"
        self._check_strategy(self.backup_strategy)
        self._check_strategy(self.default_strategy)

        self.attempted_list = set()

        self.manager_combine.set_good_pair_test(self.untried)
        self.manager_random.set_good_pair_test(self.untried)
        self.manager_untrap.set_good_pair_test(self.untried)
        self.manager_gmin.set_good_pair_test(self.untried)

    def _already_tried(self, min1, min2):

        if (min1, min2) in self.attempted_list:
            return True
        elif (min2, min1) in self.attempted_list:
            return True
        else:
            return False

    def untried(self, min1, min2):
        """return true if this minima pair have not been tried yet"""
        return not self._already_tried(min1, min2)

    def _register_pair(self, min1, min2):
        self.attempted_list.add((min1, min2))
        self.attempted_list.add((min2, min1))

    def _check_strategy(self, strategy):
        if strategy not in self.possible_strategies:
            raise Exception("strategy must be from %s" % (str(self.possible_strategies)))


    def get_connect_job(self, strategy=None):
        """return the next connect job according to the chosen strategy"""
        if strategy is None:
            strategy = self.default_strategy

        self._check_strategy(strategy)

        if strategy == "untrap":
            min1, min2 = self.manager_untrap.get_connect_job()
            if min1 is None or min2 is None:
                if self.verbosity > 0:
                    print("couldn't find any minima to untrap.  Doing", self.backup_strategy, "strategy instead")
                strategy = self.backup_strategy
            else:
                if self.verbosity > 0:
                    print("sending an untrap connect job", min1.id(), min2.id())

        if strategy == "combine":
            min1, min2 = self.manager_combine.get_connect_job()
            if min1 is None or min2 is None:
                if self.verbosity > 0:
                    print("couldn't find any minima clusters to combine.  Doing", self.backup_strategy, "strategy instead")
                strategy = self.backup_strategy
            else:
                if self.verbosity > 0:
                    print("sending a connect job to combine two disconnected clusters", min1.id(), min2.id())

        if strategy == "gmin":
            min1, min2 = self.manager_gmin.get_connect_job()
            if min1 is None or min2 is None:
                if self.verbosity > 0:
                    print("couldn't find any minima not connected to the global minimum.  Doing", self.backup_strategy, "strategy instead")
                strategy = self.backup_strategy
            else:
                if self.verbosity > 0:
                    print("sending a connect job to connect all minima with the global minimum", min1.id(), min2.id())

        if strategy == "random":
            min1, min2 = self.manager_random.get_connect_job()
            if min1 is None or min2 is None:
                raise self.NoMoreConnectionsError(
                    "couldn't find any random minima pair to connect.  Have we tried all pairs?")
            if self.verbosity > 0:
                print("sending a random connect job", min1.id(), min2.id())

        self._register_pair(min1, min2)
        return min1, min2

