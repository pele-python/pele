from pygmin.concurrent import RandomConnectWorker

nruns = 5
manager_name = "ljconnect_example"
hostname = "localhost"
port = 11567

print "I will do %d connect runs before exiting"%nruns

uri="PYRO:%s@%s:%d"%(manager_name,hostname,port)

worker = RandomConnectWorker(uri)
worker.run()