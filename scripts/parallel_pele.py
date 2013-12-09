"""
A script to run pele basinhopping and connect jobs in parallel.  This script
must know about your system, so you must define a function

    `create_system()`

which returns a pele system class object.  Save this function in a file and
pass the file name to this script.  The function will be imported using a line
analogous to `from my_system import create_system`.  We will assume
`create_system` is in file "my_system.py".  To start the server:

    python parallel_pele.py my_system.py --server
 
This will output the uri which the workers will need to connect to the server.
Now that the server is running we can submit basinhopping jobs:

    python parallel_pele.py my_system.py --basinhopping --uri=`<my_copied_uri>`

or connect jobs:

    python parallel_pele.py my_system.py --connect --uri=`<my_copied_uri>`
"""

import argparse
import Pyro4

from pele.systems import LJCluster
from pele.concurrent import ConnectServer, ConnectWorker, BasinhoppingWorker

# note: the default serializer that Pyro4 uses does not know how to 
# serialize numpy arrays so we use pickle instead.
# pickle introduces a security risk, but for this example we
# are listening only on localhost, so it is probably OK
Pyro4.config.SERIALIZER = 'pickle' 
Pyro4.config.SERIALIZERS_ACCEPTED.add('pickle')


#def create_system():
#    natoms = 38
#    system = LJCluster(natoms)
#    return system

def start_server(create_system, dbname, server_name=None, host=None, port=None):
    print "setting up system class with database", dbname
    system = create_system()
    db = system.create_database(dbname)
    
    connect_manager=ConnectServer(system, db, server_name=server_name,
                                         host=host, port=port)

    if db.number_of_minima() == 0:
        print "there are no minima in the database.  Start a basinhopping run to generate minima"
    else:
        print "Working on %d minima" % db.number_of_minima()
    
    print "to start searching for minima:"
    print "    python start_basinhopping_worker.py"
    print "to start connecting minima:"
    print "    python start_connect_worker.py"
            
    
    connect_manager.run()

def start_basinhopping(create_system, uri, niter=1000):
    system = create_system()
    worker = BasinhoppingWorker(uri, system=system)
    worker.run(nsteps=niter)

def start_connect(create_system, uri, strategy="random"):
    system = create_system()
    worker = ConnectWorker(uri, system=system, strategy=strategy)
    worker.run()



def main():
    parser = argparse.ArgumentParser(
         description=__doc__,
#         "Tool to run basinhopping and connect jobs in parallel.\n"
#         +"To use this you must first start a server."
         formatter_class=argparse.RawTextHelpFormatter
         )

    parser.add_argument("system_file", type=str, 
                        help="file containing the function `create_system()`.  "
                        "This function will be imported from this file\n"
                        "(e.g. `from <system_file> import create_system`)")
    parser.add_argument("--server", action="store_true", 
                        help="start the server")
    parser.add_argument("--database", type=str, default="db.sqlite",
                        help="database file name")
    parser.add_argument("--basinhopping", action="store_true", 
                        help="start a basinhopping run")
    parser.add_argument("--connect", action="store_true",
                        help="start a connect run")
    parser.add_argument("--uri", type=str, default=None,
                        help="uri of the server")
    parser.add_argument("--host", type=str, default="localhost",
                        help="host at which to listen for incoming messages")
    parser.add_argument("--port", type=int, default=11568,
                        help="port at which to listen for incoming messages")
    parser.add_argument("--server-name", type=str, default="pele_server",
                        help="name for the pele server")
    parser.add_argument("--strategy", type=str, default="random",
                        help="strategy to use when choosing which minima to connect")
    parser.add_argument("--niter", type=int, default=5000,
                        help="number of basinhopping iterations")
    args = parser.parse_args()
    
    if not (args.server or args.basinhopping or args.connect):
        raise ValueError("must specify either server, basinhopping, or connect")
    
    if not args.server and args.uri is None:
        raise ValueError("you must pass the URI of the pele job server")

    if args.system_file.endswith(".py"):
        args.system_file = args.system_file[:-3]
    mod = __import__(args.system_file, globals(), locals(), ["create_system"], -1)
    create_system = mod.create_system

    if args.server:
        start_server(create_system, args.database, host=args.host, port=args.port, server_name=args.server_name)
        return
    
    if args.basinhopping:
        start_basinhopping(create_system, args.uri, niter=args.niter)
        return
    
    if args.connect:
        start_connect(create_system, args.uri, strategy=args.strategy)
        return

  
if __name__ == "__main__":
    main()
