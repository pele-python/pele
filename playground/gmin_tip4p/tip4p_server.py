from random import choice
from tip4p_system import TIP4PSystem
from pele.storage import Minimum
from pele.concurrent import RandomConnectServer

server_name = "water8_connect"
hostname="localhost"
port=11569

system=TIP4PSystem()
db = system.create_database("tip4p_8.sqlite")

print "Working on %d minima"%db.session.query(Minimum).count()

print "Creating new connect manager"
connect_manager=RandomConnectServer(system, db, server_name=server_name,
                                     host=hostname, port=port)

connect_manager.run()

