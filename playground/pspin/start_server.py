"""
e.g. run python start_server.py 5 10 --server-name example --host nemesis.ch.private.cam.ac.uk
"""

import Pyro4
import argparse
import socket
import sys
import socket
import os
from pele.storage import Database
from pspin_spherical_system import  MeanFieldPSpinSphericalSystem
from pele.concurrent import ConnectServer

# note: the default serializer that Pyro4 uses does not know how to 
# serialize numpy arrays so we use pickle instead.
# pickle introduces a security risk, but for this example we
# are listening only on localhost, so it is probably OK

Pyro4.config.SERIALIZER = 'pickle'
Pyro4.config.SERIALIZERS_ACCEPTED.add('pickle')
Pyro4.config.SERVERTYPE= 'multiplex'
sys.excepthook = Pyro4.util.excepthook

def pick_unused_port():
    """
    pick an unused port number

    Returns
    -------
    port: int
        free port number
    Notes
    -----
    this does not guarantee that the port is actually free, in the short time window
    between when the port is identified as free and it is actually bound to, the port
    can be taken by another process
    """
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind(('localhost', 0))
    addr, port = s.getsockname()
    s.close()
    return port

def get_server_uri(nspins, p):
    with open('server_uri_pspin_spherical_p{}_N{}.uri'.format(p,nspins)) as f:
        uri = [line for line in f][0]
    assert uri[:5] == "PYRO:"
    return uri

def write_server_uri(server_name, hostname, port, nspins, p):
    uri = "PYRO:%s@%s:%d" % (server_name, hostname, port)
    with open('server_uri_pspin_spherical_p{}_N{}.uri'.format(p,nspins),'w') as out_server_uri:
        out_server_uri.write(uri)
    return uri

def create_system(nspins, p, interactions):
    system = MeanFieldPSpinSphericalSystem(nspins, p=p, interactions=interactions)
    return system

def _get_database_params(dbname):
    db = Database(dbname, createdb=False)
    interactions = db.get_property("interactions").value()
    db_nspins = db.get_property("nspins").value()
    db_p = db.get_property("p").value()
    params = (db_nspins, db_p, interactions)
    return db, params

def get_database_params_worker(nspins, p):
    db, (db_nspins, db_p, interactions) = _get_database_params("pspin_spherical_p{}_N{}.sqlite".format(p,nspins))
    #close this SQLAlchemy session
    db.session.close()
    #check that parameters match
    assert db_nspins == nspins
    assert db_p == p
    return interactions

def get_database_params_server(nspins, p):
    db, (db_nspins, db_p, interactions) = _get_database_params("pspin_spherical_p{}_N{}.sqlite".format(p,nspins))
    #check that parameters match
    assert db_nspins == nspins
    assert db_p == p
    return db, interactions

def main():
    parser = argparse.ArgumentParser(description="dispatcher queue")
    parser.add_argument("p", type=int, help="p-spin")
    parser.add_argument("nspins", type=int, help="number of spins")
    parser.add_argument("--server-name", type=str, help="name of the dispatcher",default=None)
    parser.add_argument("--host", type=str, help="address of the host (node on which the worker is started)",default=None)
    parser.add_argument("--port", type=int, help="port number on which the worker is started)",default=0)
    args = parser.parse_args()

    #set-up database
    nspins = args.nspins
    p = args.p
    dbname = "pspin_spherical_p{}_N{}.sqlite".format(p,nspins)
    print "setting up p={} N={} with database {}".format(p, nspins, dbname)

    #deal with existing database (if calculations has to be restarted)
    try:
        db, interactions = get_database_params_server(nspins, p)
        print "Warning: database {} already exists, using the already existing database".format(dbname)
    except IOError:
        db = None
        interactions = None

    system = create_system(nspins, p, interactions)
    if db is None:
        db = system.create_database(dbname)

    #start connect manager
    server_name = args.server_name
    if args.host is None:
        hostname = socket.gethostname()
        host = Pyro4.socketutil.getIpAddress(hostname, workaround127=True)
    else:
        host = args.host
    if args.port == 0:
        port = pick_unused_port()
    else:
        port = args.port

    connect_manager = ConnectServer(system, db, server_name=server_name,
                                    host=host, port=port)

    print "printing server uri..."
    uri = write_server_uri(server_name, host, port, nspins, p)
    print "done"

    if db.number_of_minima() == 0:
        print "there are no minima in the database.  Start a basinhopping run to generate minima"
    else:
        print "Working on %d minima" % db.number_of_minima()

    print "to start searching for minima:"
    print "    python start_basinhopping_worker.py"
    print "to start connecting minima:"
    print "    python start_connect_worker.py"

    connect_manager.run()

if __name__ == "__main__":
    main()
