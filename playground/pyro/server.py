from random import choice
from pygmin.systems import LJCluster
from pygmin.storage import Minimum
from pygmin.concurrent import RandomConnectServer

server_name = "ljconnect_example"
hostname="localhost"
port=11567

print "setting up LJ38 with database lj38.sqlite"
system = LJCluster(38)
db = system.create_database("lj38.sqlite")

if db.session.query(Minimum).count() < 10:
    print "The database is empty, run basinhopping to get some minima"
    bh = system.get_basinhopping(database=db)
    bh.run(100)

print "Working on %d minima"%db.session.query(Minimum).count()

print "Running basinhopping to generate initial database"
        
print "Creating new connect manager"
connect_manager=RandomConnectServer(system, db, server_name=server_name,
                                     host=hostname, port=port)

connect_manager.run()