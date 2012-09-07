from pygmin.storage.database import Database, Minimum
from pygmin.mindist.minpermdist_stochastic import minPermDistStochastic as minpermdist
from random import choice
import networkx as nx

from  pygmin.NEB.graph import Graph
import time
import pickle

distance = 1
npairs = 200

db = Database(db = "storage.sqlite")
    

possible_starts = db.session.query(Minimum).order_by(Minimum.energy).limit(100).all()
print len(possible_starts)

t0 = time.time()

print "generating graph"
t0 = time.time()
graph = Graph(db)
print "%.1f seconds"%(time.time() - t0)



for distance in xrange(1,10):
    print "processing distance", distance
    pairs=set()
    outfile = "coords.%d.dat"%(distance)
    while True:
        start =  choice(possible_starts)
        end = start
        tries=0
        while nx.bidirectional_dijkstra(graph.graph, start, end)[0] < distance:
            end = choice(graph.graph.neighbors(end))
            tries=+1
            if tries==1000: break
        if tries==1000:
            print "failed to move %d steps"%(distance)
            continue
        if start._id > end._id:
            start, end = end, start
        pair = (start, end) 
        if pair not in pairs:
            pairs.add(pair)
            print "pair %d/%d (%d):"%(len(pairs), npairs, distance), start._id, end._id
        if len(pairs) == npairs:
            break
        
    bench_pairs=[]
    for n,pair in zip(xrange(len(pairs)), pairs):
        dist, X1, X2 = minpermdist( pair[0].coords, pair[1].coords, niter = 100 )
        print "%d/%d (%d): distance for"%(n+1, npairs, distance), pair[0]._id, pair[1]._id, dist
        bench_pairs.append((X1.copy(), X2.copy()))
    
    pickle.dump(bench_pairs, open(outfile, "w"))