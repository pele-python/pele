from pygmin.storage.database import Database,Minimum
from sqlalchemy.sql.expression import func
from random import choice
import networkx as nx

from  pygmin.NEB.graph import Graph
import time

distance = 1

db = Database(db = "storage.sqlite")
t0 = time.time()
    
print "generating graph"
t0 = time.time()
graph = Graph(db)
print "%.1f seconds"%(time.time() - t0)

start =  choice(graph.graph.nodes())
current = start
while nx.bidirectional_dijkstra(graph.graph, start, current)[0] < distance:
    current = choice(graph.graph.neighbors(current))
    