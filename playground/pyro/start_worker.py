from pygmin.concurrent import ConnectWorker
from pygmin.systems import LJCluster

from server import create_system, port, hostname, server_name

def main():

    nruns = 5
    print "I will do %d connect runs before exiting" % nruns
    
    system = create_system()
    
    uri="PYRO:%s@%s:%d" % (server_name, hostname,port)
    worker = ConnectWorker(uri, system=system, strategy="random")
    worker.run()
    
if __name__ == "__main__":
    main()