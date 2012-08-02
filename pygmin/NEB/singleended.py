'''
Single ended searches

@author: ruehle
'''

import numpy as np
from dimer import DimerSearch
from pygmin.optimize import quench
import tstools
from pygmin.utils import zeroev

def _uphill_search(x0, search, push, push_minrms):
    ev = search.tau
        
    x1 = x0.copy()
    while True:
        x1+=push*ev
        e,g = search.potential.getEnergyGradient(x1)
        rms = np.linalg.norm(g)/np.sqrt(len(g))
        if (rms > push_minrms):
            break
    
    return quench.fire(x0, search.getEnergyGradient, tol=1e-4, maxstep=0.01)
        
def find_escape_pathes(minimum, potential, graph, ntries=1, push=1.e-2, push_minrms=1.e-2):
    print "Single ended search for minimum", minimum._id, minimum.energy
    
    search = DimerSearch(minimum.coords, potential, zeroEigenVecs=zeroev.for_cluster)
   
    for i in xrange(ntries):
        
        x_ts, energy_ts, rms, tmp = _uphill_search(minimum.coords, search, push, push_minrms)
        ret1, ret2 = tstools.minima_from_ts(potential.getEnergyGradient, x_ts, displace=1e-2)
        
        min1 = graph.addMinimum(ret1[1], ret2[0])
        min2 = graph.addMinimum(ret1[1], ret2[0])
        
        if(not min1 is minimum and not min2 is minimum):
            print "Warning in single ended search: did not find initial minimum during quench from transition state"
         
        if(min1 is min2):
            print "Warning in single ended search: downhill search from transistion state ended in same minimum"
            
        ts = graph.addTransitionState(energy_ts, x_ts, min1, min2)
        print "found transition state: ", min1._id, min2._id, min1.energy,ts.energy,min2.energy
        
        search.findNextTS()
        
        
if __name__ == "__main__":
    from graph import Graph
    from connect_min import getSetOfMinLJ
     
    natoms = 8
    pot, saveit = getSetOfMinLJ(natoms)
    graph = Graph(saveit)
 
    minima = saveit.minima()
    find_escape_pathes(minima[0], pot, graph, ntries=20)
    
    #print graph
    #for node in graph.graph.nodes():
    #    print node, node.energy
    for ts in graph.storage.transition_states():
        print ts.minimum1._id,ts.minimum2._id, "E", ts.minimum1.energy, ts.minimum2.energy, ts.minimum2.energy
        
    