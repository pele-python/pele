from random import choice
from pygmin.systems import LJCluster
from pygmin.storage import Minimum
from pygmin.concurrent import ConnectServer

server_name = "ljconnect_example"
hostname="localhost"
port=11568

def create_system():
    natoms = 38
    system = LJCluster(natoms)
    return system

def main():
    
    print "setting up LJ38 with database lj38.sqlite"
    system = create_system()
    db = system.create_database("lj38.sqlite")
    
    if db.number_of_minima() < 10:
        print "Running basinhopping to generate initial database"
        print "The database is empty, run basinhopping to get some minima"
        bh = system.get_basinhopping(database=db)
        bh.run(100)
    
    print "Working on %d minima" % db.number_of_minima()
            
    print "Creating new connect manager"
    connect_manager=ConnectServer(system, db, server_name=server_name,
                                         host=hostname, port=port)
    
    connect_manager.run()
    
if __name__ == "__main__":
    main()