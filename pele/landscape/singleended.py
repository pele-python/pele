"""
Single ended searches

@author: ruehle
"""
from __future__ import print_function
from __future__ import absolute_import

__all__ = ["find_escape_paths"]

import numpy as np

from pele.transition_states import DimerSearch, minima_from_ts, zeroEV_cluster
from pele.optimize import fire


def _uphill_search(x0, search, push, push_minrms):
    ev = search.tau
    evecs = search.get_eigenvecs(x0)
    # print len(evecs)
    x1 = x0.copy()
    while True:
        x1 += push * ev
        g = search.getOrthogonalGradient(x1, evecs)
        # print np.linalg.norm(g)
        rms = np.linalg.norm(g) / np.sqrt(len(g))
        # print "rms",rms
        if rms > push_minrms:
            # print "rms final",rms
            break
    # x1 = x0 + np.random.random(x0.shape)*0.1
    # search.tau_done=[]
    # search.x0 = x1
    # search.findNextTS()
    return fire(x1, search.getEnergyGradient, tol=1e-6)


def find_escape_paths(minimum, potential, graph, ntries=1, push=1.e-2, push_minrms=1.e-2):
    raise Exception(
        "js850> this function doesn't work anymore since changing graph.addMinimum and addTransitionState.  It needs to be overhauled")
    print("Single ended search for minimum", minimum.id(), minimum.energy)

    search = DimerSearch(minimum.coords, potential, zeroEigenVecs=zeroEV_cluster)

    for i in range(ntries):

        x_ts, energy_ts, rms, tmp = _uphill_search(minimum.coords, search, push, push_minrms)
        ret1, ret2 = minima_from_ts(potential, x_ts, displace=1e-2)

        min1 = graph.addMinimum(ret1[1], ret1[0])
        min2 = graph.addMinimum(ret2[1], ret2[0])

        if not min1 is minimum and not min2 is minimum:
            print("Warning in single ended search: did not find initial minimum during quench from transition state")

        if min1 is min2:
            print("Warning in single ended search: downhill search from transistion state ended in same minimum")

        ts = graph.addTransitionState(energy_ts, x_ts, min1, min2)
        print("found transition state: ", min1.id(), min2.id(), min1.energy, ts.energy, min2.energy)

        search.findNextTS()


if __name__ == "__main__":
    from pele.landscape import TSGraph
    from .connect_min import getSetOfMinLJ

    natoms = 8

    pot, saveit = getSetOfMinLJ(natoms)
    graph = TSGraph(saveit)

    minima = saveit.minima()
    find_escape_paths(minima[0], pot, graph, ntries=20)

    # print graph
    # for node in graph.graph.nodes():
    # print node, node.energy
    for ts in graph.storage.transition_states():
        print(ts.minimum1.id(), ts.minimum2.id(), "E", ts.minimum1.energy, ts.minimum2.energy, ts.minimum2.energy)
        
    

