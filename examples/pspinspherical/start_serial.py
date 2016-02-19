from __future__ import division
import argparse
from pele.storage import Database
from pele.systems import MeanFieldPSpinSphericalSystem


def create_system(nspins, p, interactions):
    system = MeanFieldPSpinSphericalSystem(nspins, p=p, interactions=interactions)
    return system


def get_database_params(dbname, nspins, p):
    db = Database(dbname, createdb=False)
    interactions = db.get_property("interactions").value()
    db_nspins = db.get_property("nspins").value()
    db_p = db.get_property("p").value()
    # check that parameters match
    assert db_nspins == nspins
    assert db_p == p
    return db, interactions


def make_disconnectivity_graph(system, database, fname='dg.pdf', **kwargs):
    import matplotlib.pyplot as plt
    from pele.utils.disconnectivity_graph import DisconnectivityGraph, database2graph
    graph = database2graph(database)
    dg = DisconnectivityGraph(graph, **kwargs)
    dg.calculate()
    dg.plot(linewidth=1.5)
    plt.savefig(fname)


def run_double_ended_connect(system, database, strategy='random'):
    # connect the all minima to the lowest minimum
    from pele.landscape import ConnectManager
    manager = ConnectManager(database, strategy=strategy)
    for i in xrange(database.number_of_minima()-1):
        min1, min2 = manager.get_connect_job()
        connect = system.get_double_ended_connect(min1, min2, database)
        connect.connect()


def main():
    parser = argparse.ArgumentParser(description="do nested sampling on a 3 layer elu neural network")
    # pspin variables
    parser.add_argument("p", type=int, help="p-spin")
    parser.add_argument("nspins", type=int, help="number of spins")
    # operations to perform
    parser.add_argument("--bh", type=int, help="number of basin hopping steps to perform", default=0)
    parser.add_argument("--connect", action="store_true", help="run all")
    parser.add_argument("--connect-method", type=str, help="method used to connect", default='random')

    args = parser.parse_args()
    print args

    # pspin parameters
    p, nspins = args.p, args.nspins

    # operations
    bh_niter = args.bh
    connect, connect_method = args.connect, args.connect_method
    dbname = "pspin_spherical_p{}_N{}.sqlite".format(p,nspins)
    try:
        db, interactions = get_database_params(dbname, nspins, p)
        print "Warning: database {} already exists, using the already existing database".format(dbname)
    except IOError:
        db = None
        interactions = None

    system = create_system(nspins, p, interactions)

    if db is None:
        db = system.create_database(dbname)

    # now actually run the computattion
    fname = None
    if bh_niter > 0:
        bh = system.get_basinhopping(database=db, outstream=None)
        bh.run(bh_niter)
    if connect:
        fname = "pspin_spherical_p{}_N{}.dg.pdf".format(p,nspins)
        run_double_ended_connect(system, db, strategy=connect_method)

    if fname is not None:
        make_disconnectivity_graph(system, db, fname=fname)


if __name__ == "__main__":
    main()