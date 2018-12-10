"""Wrapper to represent a storage class as a graph"""
from __future__ import print_function
import networkx as nx

__all__ = ["TSGraph", "Graph", "database2graph"]


def database2graph(db, Emax=None):
    """
    make a networkx graph from a database

    Parameters
    ----------
    db : pele Database
    Emax : float optional
        including only transition states with energy < Emax
    """
    from pele.storage.database import Minimum, TransitionState

    g = nx.Graph()
    # js850> It's not strictly necessary to add the minima explicitly here,
    # but for some reason it is much faster if you do (factor of 2).  Even 
    # if this means there are many more minima in the graph.  I'm not sure 
    # why this is.  This step is already often the bottleneck of the d-graph 
    # calculation.
    if Emax is not None:
        minima = db.session.query(Minimum).filter(Minimum.energy <= Emax)
    else:
        minima = db.session.query(Minimum)
    g.add_nodes_from(minima)
    # if we order by energy first and add the transition states with the largest energy first
    # then we will take the smallest energy transition state in the case of duplicates
    if Emax is not None:
        ts = db.session.query(TransitionState).filter(TransitionState.energy <= Emax) \
            .order_by(-TransitionState.energy)
    else:
        ts = db.session.query(TransitionState).order_by(-TransitionState.energy)
    for t in ts:
        g.add_edge(t.minimum1, t.minimum2, ts=t)
    return g


class _ConnectedComponents(nx.utils.UnionFind):
    """
    a class to build and maintain the connected components
    
    This allows connections to be determined much more rapidly because
    the breadth first search algorithm is slow
    """

    def are_connected(self, m1, m2):
        return self[m1] == self[m2]


class TSGraph(object):
    """
    Wrapper to represent a database object as a graph

    This class is primarily used by DoubleEndedConnect and has a number
    of utility functions for that purpose.  If you want access to the
    networkx Graph object, it is stored as self.graph

    Parameters
    ----------
    database :
        the database object to represent
    minima : list of minima, optional
        if none, include all minima and transition states from the database,
        else include only the minima in the list
    no_edges : bool, optional
        If true include no transition states in the graph.
        This is used, for example, to try to find a new connection between
        two minima.

    See Also
    --------
    DoubleEndedConnect

    Examples
    --------
    a graph can be easily constructed from a database::

    >>> graph = TSGraph(database)

    the networkx graph is accessed directly by

    >>> networkx_graph = graph.graph
    """

    def __init__(self, database, minima=None, no_edges=False):
        self.graph = nx.Graph()
        self.storage = database
        self.connected_components = _ConnectedComponents()
        self.minima = minima
        self.no_edges = no_edges
        self.refresh()

    def _build_all(self):
        """
        add all minima and all transition states to the graph
        """
        for m in self.storage.minima():
            self.graph.add_node(m)
        if not self.no_edges:
            for ts in self.storage.transition_states():
                self.graph.add_edge(ts.minimum1, ts.minimum2, ts=ts)
                self.connected_components.union(ts.minimum1, ts.minimum2)

    def _build_from_list(self, minima):
        """
        add only those minima from the list `minima` to the graph.
        Don't add any transition states
        """
        minima = set(minima)
        for m in minima:
            self.graph.add_node(m)
        if not self.no_edges:
            for ts in self.storage.transition_states():
                m1, m2 = ts.minimum1, ts.minimum2
                if m1 in minima:
                    if m2 in minima:
                        self.graph.add_edge(ts.minimum1, ts.minimum2, ts=ts)
                        self.connected_components.union(ts.minimum1, ts.minimum2)

    def refresh(self):
        if self.minima is None:
            self._build_all()
        else:
            self._build_from_list(self.minima)

    def addMinimum(self, minimum):
        """
        add a minimum to the database and graph
        """
        self.graph.add_node(minimum)
        return minimum

    def addTransitionState(self, ts):
        self.graph.add_edge(ts.minimum1, ts.minimum2, ts=ts)
        self.connected_components.union(ts.minimum1, ts.minimum2)
        return ts

    def areConnected(self, min1, min2):
        return self.connected_components.are_connected(min1, min2)

    def getPath(self, min1, min2):
        try:
            return nx.bidirectional_dijkstra(self.graph, min1, min2)
        except nx.NetworkXNoPath:
            return None

    def mergeMinima(self, min1, min2):
        """
        delete minima2.  all transition states pointing to min2 should
        now point to min1
        """
        self.connected_components.union(min1, min2)
        # make the edges of min2 now point to min1
        for v, data in self.graph[min2].items():
            if v == min1: continue
            if v == min2:
                print("warning: minimum", min2.id(), "is connected to itself")
                continue
            # the new edge will be (min1, v).  Add it if it doesn't already exist
            if not self.graph.has_edge(min1, v):
                self.graph.add_edge(min1, v, **data)
        self.graph.remove_node(min2)


class Graph(TSGraph):
    """this is included for backwards compatibility"""
    pass

