from pygmin.concurrent import ConnectWorker
from pygmin.systems import LJCluster

nruns = 5
manager_name = "ljconnect_example"
hostname = "localhost"
port = 11568

print "I will do %d connect runs before exiting"%nruns

uri="PYRO:%s@%s:%d"%(manager_name,hostname,port)

system = LJCluster(38)


worker = ConnectWorker(uri, system=system, strategy="combine")
worker.run()