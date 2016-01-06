import argparse
from pele.concurrent import BasinhoppingWorker
from start_server import create_system, get_server_uri
from pele.storage import Database

def main():
    parser = argparse.ArgumentParser(description="connect worker queue")
    parser.add_argument("p", type=int, help="p-spin")
    parser.add_argument("nspins", type=int, help="number of spins")
    parser.add_argument("--nsteps", type=int, help="number of basin hopping steps", default=1000)
    args = parser.parse_args()

    nspins = args.nspins
    p = args.p

    #get interactions from database
    dbname = "pspin_spherical_p{}_N{}.sqlite".format(p,nspins)
    db = Database(dbname, createdb=False)
    interactions = db.get_property("interactions").value()
    db_nspins = db.get_property("nspins").value()
    db_p = db.get_property("p").value()
    assert db_nspins == nspins
    assert db_p == p
    #close this SQLAlchemy session
    db.session.close()

    system = create_system(nspins, p, interactions)

    uri = get_server_uri()
    worker = BasinhoppingWorker(uri, system=system)
    worker.run(args.nsteps)

if __name__ == "__main__":
    main()
