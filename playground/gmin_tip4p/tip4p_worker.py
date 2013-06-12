from tip4p_system import TIP4PSystem
from pele.concurrent import RandomConnectWorker

manager_name = "water8_connect"
hostname = "localhost"
port = 11569

uri="PYRO:%s@%s:%d"%(manager_name,hostname,port)

system=TIP4PSystem()
worker = RandomConnectWorker(uri, system=system)
worker.run()

