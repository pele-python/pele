import numpy as np
from pygmin.NEB.NEB import NEB
from pygmin.optimize.transition_state.transition_state_refinement import findTransitionState
import tstools
from pygmin.storage.savenlowest import Minimum 
from min_ts_graph import TransitionState
import networkx as nx

def connectMinima(min1, min2, pot, graph, mindist):
    """
    create a connection between two minima composed of discrete sets of minima and transition states
    """
    min1id = graph.addMin(min1)
    min2id = graph.addMin(min2)
    return connectMinimaID(min1id, min2id, pot, graph, mindist)

def connectMinimaID(min1id, min2id, pot, graph, mindist):
    if graph.areConnected(min1id, min2id):
        return
    print "warning connectMinima is still in beta"

    print ""
    print "starting connect run to try to connect minima", min1id, min2id
    print "minimizing the distance between the two minima"
    min1 = graph.node2min[min1id]
    min2 = graph.node2min[min2id]
    dist, newcoords1, newcoords2 = mindist(min1.coords, min2.coords) 
    print "runing NEB"
    neb = NEB(newcoords1, newcoords2, pot)
    neb.optimize()
    neb.MakeAllMaximaClimbing()
    neb.optimize()

    minimalist = [min1id]
    nclimbing = np.sum( [neb.isclimbing[i] for i in range(neb.nimages) ])
    print "from NEB search found", nclimbing, "climbing images"
    
    for i in range(neb.nimages):
        if neb.isclimbing[i]:
            coords = neb.coords[i,:]
            print "refining transition state from NEB climbing image"
            ret = findTransitionState(coords, pot)
            coords, eigval, eigvec, E, grad, rms = ret

            ts = TransitionState(coords, E, eigvec, eigval)
            print "falling off either side of transition state to find new minima"
            ret1, ret2 = tstools.minima_from_ts(pot.getEnergyGradient, coords, \
                displace=1e-1, quenchParameters={"tol":1e-7})
            
            newmin1 = Minimum(ret1[1], ret1[0])
            newmin2 = Minimum(ret2[1], ret2[0])
            id1 = graph.addMin(newmin1)
            id2 = graph.addMin(newmin2)
            if id1 == id2:
                print "warning: stepping off the transition state resulted in twice the same minima", id1
                continue
            print "adding transition state", id1, id2
            graph.addTS(id1, id2, ts)
            

            
            connected = graph.areConnected(min1id, min2id)
            print "connected yet?", connected
            if connected:
                return minimalist
            
            if id2 == min2id: continue
            if id1 == min2id: continue
            if not id1 in minimalist:
                minimalist.append(id1)
            if not id2 in minimalist:
                minimalist.append(id2)

    minimalist.append(min2id)
    
    print ""
    print "done with NEB, the minima are still not connected."
    print "will now try to fill in the gaps"
    print "the minima to connect are", minimalist
    #we should choose the next minima pair to try to connect in a more clever way.
    for i in range(len(minimalist)-1):
        id1 = minimalist[i]
        id2 = minimalist[i+1]
        print "id1 id2", id1, id2
        #note: this should be redone, it can enter an endless loop in some cases
        #  if the order goes something like this
        #  connect(1,2)
        #      connect(1,3)
        #          connect(1,2)
        # minimalist should be passed to connect
        ret = connectMinimaID(id1, id2, pot, graph, mindist)
        if graph.areConnected(min1id, min2id):
            return






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
 
    connectMinima(min1, min2, pot, graph, mindist)
    
    print graph
    for node in graph.graph.nodes():
        print node, graph.node2min[node].E
    for edge, ts in graph.tstates().items():
        node1, node2 = edge 
        E1 = graph.node2min[node1].E
        E2 = graph.node2min[node2].E
        Ets = ts.E
        print node1, node2, "E", E1, Ets, E2
        
    distances, path = nx.bidirectional_dijkstra(graph.graph, min1id, min2id)
    with open("path.out", "w") as fout:
        for i in range(len(path)-1):
            n1 = path[i]
            n2 = path[i+1]
            m1 = graph.node2min[n1]
            m2 = graph.node2min[n2]
            ts = graph.getTS(n1, n2)
            print "path", n1, "->", n2, m1.E, "/->", ts.E, "\->", m2.E
            fout.write("%f\n" % m1.E)
            fout.write("%f\n" % ts.E)
        n2 = path[-1]
        m2 = graph.node2min[n2]
        fout.write("%f\n" % m2.E)


if __name__ == "__main__":
    test()


    