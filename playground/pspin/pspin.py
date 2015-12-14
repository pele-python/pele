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
import cmath

def isclose(a, b, rel_tol=1e-9, abs_tol=0.0, method='weak'):
    """
    code imported from math.isclose python 3.5
    """
    if method not in ("asymmetric", "strong", "weak", "average"):
        raise ValueError('method must be one of: "asymmetric",'
                         ' "strong", "weak", "average"')

    if rel_tol < 0.0 or abs_tol < 0.0:
        raise ValueError('error tolerances must be non-negative')
    
    if a == b:  # short-circuit exact equality
        return True
    # use cmath so it will work with complex or float
    if cmath.isinf(a) or cmath.isinf(b):
        # This includes the case of two infinities of opposite sign, or
        # one infinity and one finite number. Two infinities of opposite sign
        # would otherwise have an infinite relative tolerance.
        return False
    diff = abs(b - a)
    if method == "asymmetric":
        return (diff <= abs(rel_tol * b)) or (diff <= abs_tol)
    elif method == "strong":
        return (((diff <= abs(rel_tol * b)) and
                 (diff <= abs(rel_tol * a))) or
                (diff <= abs_tol))
    elif method == "weak":
        return (((diff <= abs(rel_tol * b)) or
                 (diff <= abs(rel_tol * a))) or
                (diff <= abs_tol))
    elif method == "average":
        return ((diff <= abs(rel_tol * (a + b) / 2) or
                (diff <= abs_tol)))
    else:
        raise ValueError('method must be one of:'
                         ' "asymmetric", "strong", "weak", "average"')

def compare_exact(x1, x2, 
                  rel_tol=1e-9,
                  abs_tol=0.0,
                  method='weak'):
    assert x1.size == x2.size
    N = x1.size
    assert isclose(np.dot(x1,x1), N)
    assert isclose(np.dot(x2,x2), N)
    dot = np.dot(x1, x2)
    return isclose(dot, N, rel_tol=rel_tol, abs_tol=abs_tol, method=method)

def func(pot, coords):
    #print coords
    #print "start energy", pot.getEnergy(coords)
    results = lbfgs_cpp(coords, pot, nsteps=1e5, tol=1e-9, iprint=-1, maxstep=10)
    #results = modifiedfire_cpp(coords, pot, nsteps=1e5, tol=1e-5, iprint=-1)
    #print "quenched energy", results.energy
    if results.success:
        return [results.coords, results.energy, results.nfev]    

def main():
    p=3
    nspins=20
    interactions = np.ones(np.power(nspins,p))
    coords = np.ones(nspins)
    pot = MeanFieldPSpinSpherical(interactions, nspins, p, tol=1e-6)
    e = pot.getEnergy(coords)
    assert e + comb(nspins,p)/np.power(nspins,(p-1)/2) < 1e-10
    print "passed"
    
    #interactions = np.random.normal(0, np.sqrt(factorial(p)), [nspins for i in xrange(p)])
    assert p==3, "the interaction matrix setup at the moment requires that p==3"
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
    
    for _ in xrange(1000):
        coords = np.random.normal(0, 1, nspins)
        coords /= (np.linalg.norm(coords)/np.sqrt(nspins))
        coords_list.append(coords)
    
    if False:
        
        start = time.time()
        
        out = Parallel(n_jobs=max(1,4))(delayed(func)(pot, x) for x in coords_list)
        out = np.array(out)
        
        done = time.time()
        elapsed = done - start
        print 'elapsed time: ',elapsed
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(out[:,1]/nspins)
        ax.set_xlabel('E/N')
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        ax1.hist(out[:,2])
        ax1.set_xlabel('nfev')
        fig.savefig('energy_histogram_n{}.pdf'.format(nspins))
        fig1.savefig('nfev_histogram_n{}.pdf'.format(nspins))
        plt.show()
    
    if True:
        #check how many different minima
        
        start = time.time()
        
        out = Parallel(n_jobs=max(1,4))(delayed(func)(pot, x) for x in coords_list)
        out = np.array(out)
        done = time.time()
        elapsed = done - start
        print 'elapsed time: ',elapsed
        
        #print out[:,0]
        
        uniquex = []
        for x1 in out[:,0]:
            unique = True
            for x2 in uniquex:
                if compare_exact(x1, x2, rel_tol=1e-7):
                    unique = False
                    break
            if unique:
                uniquex.append(x1)
                 
        print "distinct minima", len(uniquex)
        
    
    if False:
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
    
