import matplotlib.pyplot as plt

import pygmin.utils.disconnectivity_graph as dg
from pygmin.storage import Database
from pygmin.landscape import Graph
from pygmin.systems import LJCluster

def get_database(natoms=13, nconn=5):
    """create a database for a lennard jones system
    
    fill it with minima from a basinhopping run, then connect
    some of those minima using DoubleEndedConnect
    """
    ljsys = LJCluster(natoms)
    db = ljsys.create_database()
    bh = ljsys.get_basinhopping(database=db, outstream=None)
    while (len(db.minima()) < nconn+1):
        bh.run(100)
    
    
    minima = list(db.minima())
    m1 = minima[0]
    for m2 in minima[1:nconn+1]:
        connect = ljsys.get_double_ended_connect(m1, m2, db)
        connect.connect()
    
    return db



def make_graph(database):
    #make a graph from the database
    graphwrapper = Graph(db)
    
    #turn the graph into a disconnectivity graph
    mydg = dg.DisconnectivityGraph(graphwrapper.graph, 
                                   nlevels=5,
                                   center_gmin=False,
                                   order_by_energy=True,
#                                   Emax=-169.
                                   )
    mydg.calculate()
    
    print "number of minima:", mydg.tree_graph.number_of_leaves()
    mydg.plot()
    plt.show()

if __name__ == "__main__":
    if True:
        db = get_database()
    else:
        db = Database("lj38.sqlite")
        #db = Database("database.sqlite.large")

    make_graph(db)