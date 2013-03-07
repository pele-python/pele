from random import choice
import Pyro4
from pygmin.systems import LJCluster
from pygmin.storage import Minimum

manager_name = "ljconnect_example"
hostname="localhost"
port=11567

# we need to run pyros in multiplex mode, otherwise we run into problems with 
# SQLAlchemy. This is due to the fact that a session can only be used from one
# thread. We really should fix this issue and allow for multiple sessions!
Pyro4.config.SERVERTYPE = "multiplex"


class ConnectManager(object):
    ''' Manager which decides which connect jobs to run '''
    
    def __init__(self, system, database):
        self.system = system
        self.db = database
        
    def get_connect_job(self):
        ''' get a new connect job '''
        print "new connect job assigned to worker"
        minima = self.db.minima()
        min1 = choice(minima)
        min2 = choice(minima)
        
        return min1._id, min1.coords, min2._id, min2.coords

    def get_system(self):
        ''' get the system class '''
        return self.system
    
    def add_minimum(self, E, coords):
        ''' called by worker if a new minimum is found '''
        print "a client found a minimum", E
        return self.db.addMinimum(E, coords)._id
    
    def add_ts(self, id1, id2, E, coords):
        ''' called by worker if a new transition state is found '''
        
        print "a client found a transition state", E
        min1 = self.db.session.query(Minimum).filter(Minimum._id == id1).one()
        min2 = self.db.session.query(Minimum).filter(Minimum._id == id2).one()
        
        return self.db.addTransitionState(E, coords, min1, min2)._id

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
connect_manager=ConnectManager(system, db)

print "Starting Pyros daemon"
daemon=Pyro4.Daemon(host=hostname, port=port)

# make the connect_manager available to Pyros childs
uri=daemon.register(connect_manager, objectId=manager_name)
print "The connect manager can be accessed by the following uri: ", uri 

print "Ready to accept connections"
daemon.requestLoop() 