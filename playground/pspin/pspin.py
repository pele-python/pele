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

def main():
    p=3
    nspins=50
    interactions = np.ones(np.power(nspins,p))
    coords = np.ones(nspins)
    pot = MeanFieldPSpinSpherical(interactions, nspins, p, tol=1e-6)
    e = pot.getEnergy(coords)
    assert e == -comb(nspins,p)/nspins
    
    interactions = np.random.normal(0, np.sqrt(factorial(p)), [nspins for i in xrange(p)])
    for i in xrange(nspins):
        for j in xrange(i, nspins):
            for k in xrange(j, nspins):
                interactions[k][i][j] = interactions[i][j][k]
                interactions[k][j][i] = interactions[i][j][k]
                interactions[j][k][i] = interactions[i][j][k]
                interactions[i][k][j] = interactions[i][j][k]
                interactions[j][i][k] = interactions[i][j][k]
    
    interactionsfl = interactions.flatten()
    pot = MeanFieldPSpinSpherical(interactionsfl, nspins, p, tol=1e-6)
    energies = []
    nfevs = []
    coords_list = []
      
    for _ in xrange(1):
        coords = np.random.normal(0, 1, nspins)
        coords /= (np.linalg.norm(coords)/np.sqrt(nspins))
        coords_list.append(coords)
    
    if False:
        def func(coords):
            #print coords
            print "start energy", pot.getEnergy(coords)
            results = lbfgs_cpp(coords, pot, nsteps=1e5, tol=1e-5, iprint=-1, maxstep=10)                                                                                                                                 
            #results = modifiedfire_cpp(coords, pot, nsteps=1e6, tol=1e-5, iprint=-1)
            print "quenched energy", results.energy
            if results.success:
                return [results.energy, results.nfev]
        
        out = Parallel(n_jobs=max(1,7))(delayed(func)(x) for x in coords_list)
        out = np.array(out)
        
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
        G = nx.Graph()
        G.add_nodes_from(xrange(nspins))
        
        com = combinations(range(nspins), p)
        l = 0
        for c in com:
            i, j, k = c
            if i != j and j != k and i != k:
                w = interactions[i][j][k]
                G.add_edge(i, j, weight=w)
                G.add_edge(j, k, weight=w)
                G.add_edge(k, i, weight=w)
            l += 1
        assert comb(nspins,p) == l
        
        # use one of the edge properties to control line thickness
        epos = [(u,v) for (u,v,d) in sorted(G.edges(data=True), key = lambda (a, b, dct): dct['weight'], reverse=True)]
        eneg = [(u,v) for (u,v,d) in sorted(G.edges(data=True), key = lambda (a, b, dct): dct['weight'])]
        epos = epos[:min(len(epos),min(nspins*2,50))]
        eneg = eneg[:min(len(eneg),min(nspins*2,50))]
        wei = nx.get_edge_attributes(G, 'weight')
        print [wei[(u,v)] for (u,v) in epos]
        print [wei[(u,v)] for (u,v) in eneg]
        # layout
        pos = nx.spring_layout(G, iterations=int(1e3))
        #pos = nx.circular_layout(G)
        
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111)
        plt.axis('off')
        
        # rendering
        nodesize = abs(coords_list[0])/np.amax(coords_list[0])*1e3
        nx.draw_networkx_nodes(G, pos, node_color=coords_list[0], cmap=plt.cm.get_cmap('RdYlBu'), alpha=0.7, 
                               linewidths=0, node_size=nodesize, label=range(nspins), ax=ax)
        nx.draw_networkx_edges(G, pos, edgelist=epos, width=2, alpha=0.3, edge_color='r', ax=ax)
        nx.draw_networkx_edges(G, pos, edgelist=eneg, width=2, alpha=0.3, edge_color='b', style='dashed', ax=ax)
        nx.draw_networkx_labels(G,pos, font_color='', ax=ax)
        
        plt.axes().set_aspect('equal', 'datalim')
        #plt.show()
        plt.savefig('pspin_n{}_network.pdf'.format(nspins))

if __name__ == "__main__":
    main()
    
