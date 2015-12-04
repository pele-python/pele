from __future__ import division
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.misc import comb, factorial
from pele.optimize._quench import lbfgs_cpp, modifiedfire_cpp
from pele.potentials import MeanFieldPSpinSpherical
from joblib import Parallel, delayed
import networkx as nx
from itertools import combinations
import time

def func(pot, coords):
    #print coords
    print "start energy", pot.getEnergy(coords)
    results = lbfgs_cpp(coords, pot, nsteps=1e5, tol=1e-5, iprint=-1, maxstep=10)                                                                                                                                 
    #results = modifiedfire_cpp(coords, pot, nsteps=1e6, tol=1e-5, iprint=-1)
    print "quenched energy", results.energy
    if results.success:
        return [results.energy, results.nfev]

def main():
    p=3
    nspins=250
    interactions = np.ones(np.power(nspins,p))
    coords = np.ones(nspins)
    pot = MeanFieldPSpinSpherical(interactions, nspins, p, tol=1e-6)
    e = pot.getEnergy(coords)
    assert e + comb(nspins,p)/np.power(nspins,(p-1)/2) < 1e-10
    print "passed"
    
    #interactions = np.random.normal(0, np.sqrt(factorial(p)), [nspins for i in xrange(p)])
    interactions = np.empty([nspins for i in xrange(p)])
    for i in xrange(nspins):
        for j in xrange(i, nspins):
            for k in xrange(j, nspins):
                w = np.random.normal(0, np.sqrt(factorial(p)))
                interactions[i][j][k] = w
                interactions[k][i][j] = w
                interactions[k][j][i] = w
                interactions[j][k][i] = w
                interactions[i][k][j] = w
                interactions[j][i][k] = w
    
    interactionsfl = interactions.flatten()
    pot = MeanFieldPSpinSpherical(interactionsfl, nspins, p, tol=1e-6)
    energies = []
    nfevs = []
    coords_list = []
    
    print "here1"  
    
    for _ in xrange(1):
        coords = np.random.normal(0, 1, nspins)
        coords /= (np.linalg.norm(coords)/np.sqrt(nspins))
        coords_list.append(coords)
    
    print "here2"
    
    if False:
        
        start = time.time()
        
        out = Parallel(n_jobs=max(1,7))(delayed(func)(pot, x) for x in coords_list)
        out = np.array(out)
        
        done = time.time()
        elapsed = done - start
        print 'elapsed time: ',elapsed
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(out[:,0]/nspins)
        ax.set_xlabel('E/N')
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        ax1.hist(out[:,1])
        ax1.set_xlabel('nfev')
        fig.savefig('energy_histogram_n{}.pdf'.format(nspins))
        fig1.savefig('nfev_histogram_n{}.pdf'.format(nspins))
        plt.show()
    
    if True:
        # create a graph object, add n nodes to it, and the edges
        Gm = nx.MultiGraph()
        Gm.add_nodes_from(xrange(nspins))
        
        l = 0
        for c in combinations(range(nspins), p):
            i, j, k = c
            if i != j and j != k and i != k:
                w = interactions[i][j][k]
                Gm.add_edge(i, j, weight=w)
                Gm.add_edge(j, k, weight=w)
                Gm.add_edge(k, i, weight=w)
            l += 1
        print l*3
        assert comb(nspins,p) == l
        assert l*3 == len(Gm.edges())
        
        #merge edges of multigraph
        G = nx.Graph()
        for u,v,data in Gm.edges_iter(data=True):
            w = data['weight']
            if G.has_edge(u,v):
                G[u][v]['weight'] += w
            else:
                G.add_edge(u, v, weight=w)
        
        # use one of the edge properties to control line thickness
        epos = [(u,v) for (u,v,d) in sorted(G.edges(data=True), key = lambda (a, b, dct): dct['weight'], reverse=True) if d['weight'] > 0]
        eneg = [(u,v) for (u,v,d) in sorted(G.edges(data=True), key = lambda (a, b, dct): dct['weight']) if d['weight'] <= 0]
        assert len(epos)+len(eneg) == len(G.edges())
#        print epos
#        print [(u, v, d['weight']) for (u, v, d) in sorted(G.edges(data=True), key = lambda (a, b, dct): dct['weight'], reverse=True) if d['weight'] > 0]
#        print "eneg"
#        print eneg
#        print [(u, v, d['weight']) for (u, v, d) in sorted(G.edges(data=True), key = lambda (a, b, dct): dct['weight']) if d['weight'] <= 0]
        maxlen = 50 - 50 % p
        epos = epos[0:min(len(epos),maxlen)]
        eneg = eneg[0:min(len(eneg),maxlen)]
#        print "truncated"
#        for (u, v) in epos:
#            print u, v, G[u][v]['weight']
        
        #print [wei[(u,v)] for (u,v) in eneg]
        # layout
        pos = nx.spring_layout(G, iterations=int(1e3))
        #pos = nx.circular_layout(G)
        
        fig = plt.figure(figsize=(10,20))
        ax = fig.add_subplot(211)
        
        plt.axis('off')
        ax.set_title('Random coords')
        # rendering
        nodesize = abs(coords_list[0])/np.amax(coords_list[0])*1e3
        nx.draw_networkx_nodes(G, pos, node_color=coords_list[0], cmap=plt.cm.get_cmap('RdYlBu'), alpha=0.7, 
                               linewidths=0, node_size=nodesize, label=range(nspins), ax=ax)
        nx.draw_networkx_edges(G, pos, edgelist=epos, width=2, alpha=0.3, edge_color='r', ax=ax)
        nx.draw_networkx_edges(G, pos, edgelist=eneg, width=2, alpha=0.3, edge_color='b', style='dashed', ax=ax)
        #nx.draw_networkx_labels(G, pos, font_color='k', ax=ax)
        
        ax2 = fig.add_subplot(212)
        ax2.set_title('Minimum coords')
        plt.axis('off')
        
        # rendering
        minimum = lbfgs_cpp(coords, pot, nsteps=1e5, tol=1e-5, iprint=10, maxstep=10).coords
        nodesize = abs(minimum)/np.amax(minimum)*1e3
        nx.draw_networkx_nodes(G, pos, node_color=minimum, cmap=plt.cm.get_cmap('RdYlBu'), alpha=0.7, 
                               linewidths=0, node_size=nodesize, label=range(nspins), ax=ax2)
        nx.draw_networkx_edges(G, pos, edgelist=epos, width=2, alpha=0.3, edge_color='r', ax=ax2)
        nx.draw_networkx_edges(G, pos, edgelist=eneg, width=2, alpha=0.3, edge_color='b', style='dashed', ax=ax2)
        #nx.draw_networkx_labels(G, pos, font_color='k', ax=ax2)
        
        #plt.axes().set_aspect('equal', 'datalim')
        plt.savefig('pspin_n{}_network.pdf'.format(nspins))
        #plt.show()

if __name__ == "__main__":
    main()
    
