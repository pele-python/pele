from pele.landscape import database2graph
from pele.utils.disconnectivity_graph import DisconnectivityGraph
from pele.storage import Database
from pele.systems import LJCluster


def get_database(natoms=13, nconn=5):
    """create a database for a lennard jones system
    
    fill it with minima from a basinhopping run, then connect
    some of those minima using DoubleEndedConnect
    """
    ljsys = LJCluster(natoms)
    db = ljsys.create_database()
    bh = ljsys.get_basinhopping(database=db, outstream=None)
    while (len(db.minima()) < nconn + 1):
        bh.run(100)

    minima = list(db.minima())
    m1 = minima[0]
    for m2 in minima[1:nconn + 1]:
        connect = ljsys.get_double_ended_connect(m1, m2, db)
        connect.connect()

    return db


def make_graph(database):
    # make a graph from the database
    graph = database2graph(database)

    # turn the graph into a disconnectivity graph
    dg = DisconnectivityGraph(graph,
                              nlevels=5,
                              center_gmin=False,
                              order_by_energy=True)
    dg.calculate()

    print "number of minima:", dg.tree_graph.number_of_leaves()
    dg.plot()
    dg.show()


if __name__ == "__main__":
    if True:
        db = get_database()
    else:
        db = Database("lj38.sqlite")

    make_graph(db)
