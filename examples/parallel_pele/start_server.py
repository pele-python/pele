import Pyro4

from pele.systems import LJCluster
from pele.concurrent import ConnectServer

# note: the default serializer that Pyro4 uses does not know how to 
# serialize numpy arrays so we use pickle instead.
# pickle introduces a security risk, but for this example we
# are listening only on localhost, so it is probably OK
Pyro4.config.SERIALIZER = 'pickle'
Pyro4.config.SERIALIZERS_ACCEPTED.add('pickle')

server_name = "ljconnect_example"
hostname = "esperanto"
port = 11568


def get_server_uri():
    uri = "PYRO:%s@%s:%d" % (server_name, hostname, port)
    return uri


def create_system():
    natoms = 38
    system = LJCluster(natoms)
    return system


def main():
    print "setting up LJ38 with database lj38.sqlite"
    system = create_system()
    db = system.create_database("lj38.sqlite")

    connect_manager = ConnectServer(system, db, server_name=server_name,
                                    host=hostname, port=port)

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
