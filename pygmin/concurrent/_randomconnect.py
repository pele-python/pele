from random import choice
from collections import deque

import Pyro4

from pygmin.storage import Minimum, TransitionState
from pygmin.landscape import ConnectManager


__all__ = ["RandomConnectServer", "RandomConnectWorker"]

# we need to run pyros in multiplex mode, otherwise we run into problems with 
# SQLAlchemy. This is due to the fact that a session can only be used from one
# thread. We really should fix this issue and allow for multiple sessions!
Pyro4.config.SERVERTYPE = "multiplex"
        

class RandomConnectServer(object):
    ''' 
    Manager which decides which connect jobs to run 
    
    The manager also receives minima and transition states from the workers
    and adds them to the database. 
    
    Parameters
    ----------
    system : pygmin.system.BaseSystem
        system class to process
    database : pygmin.storage.Database
        working database
    server_name : string, optional
        Unique name for clients to connect to this server on current host 
        (objid for pyros). None for random
    host : string, optional
        host to setup server. default is localhost which does not allow connections
        from remote machines
    port : integer, optional
        port to listen for connections
    '''
    
    def __init__(self, system, database, server_name=None, host=None, port=0):
        self.system = system
        self.db = database
        self.manager_name = server_name
        self.host=host
        self.port=port
        self.Emax = None
        
        self.list_len = 100
        
        self.connect_manager = ConnectManager(self.db, Emax=self.Emax)
        
        
    def set_emax(self, Emax):
        raise Exception("set_emax is not implemented yet in the new ConnectManager scheme")
        self.Emax = None

    def get_connect_job(self, strategy="random"):
        ''' get a new connect job '''
        min1, min2 = self.connect_manager.get_connect_job(strategy)
        return min1._id, min1.coords, min2._id, min2.coords

    def get_system(self):
        ''' provide system class to worker '''
        return self.system
    
    def add_minimum(self, E, coords):
        ''' called by worker if a new minimum is found '''
        print "a client found a minimum", E
        return self.db.addMinimum(E, coords)._id
    
    def add_ts(self, id1, id2, E, coords, eigenval=None, eigenvec=None):
        ''' called by worker if a new transition state is found '''
        
        print "a client found a transition state", E
        min1 = self.db.session.query(Minimum).filter(Minimum._id == id1).one()
        min2 = self.db.session.query(Minimum).filter(Minimum._id == id2).one()
        
        return self.db.addTransitionState(E, coords, min1, min2, eigenval=eigenval, eigenvec=eigenvec)._id

    def run(self):
        ''' start the server and listen for incoming connections '''
        print "Starting Pyros daemon"
        daemon=Pyro4.Daemon(host=self.host, port=self.port)
        # make the connect_manager available to Pyros childs
        uri=daemon.register(self, objectId=self.manager_name)
        print "The connect manager can be accessed by the following uri: ", uri 
        
        print "Ready to accept connections"
        daemon.requestLoop() 
        
class RandomConnectWorker(object):
    ''' 
    worker class to execute connect runs.
    
    The worker will return all minima and transiton states found to the server
    
    Parameters
    ----------
    uri : string
        uri for job server
    system : BaseSystem, optional
        if no system class is specified, the worker obtains the system
        class from the Manager by get_system. This only works for pickleable
        systems classes. If this is not the case, the system class can be
        created on the client side and passed as a parameter.
    '''
    
    def __init__(self,uri, system=None, strategy="random"):
        print "connecting to",uri
        self.connect_manager=Pyro4.Proxy(uri)
        if system is None:
            system = self.connect_manager.get_system()
        self.system = system
        
        self.strategy = strategy
        
    def run(self, nruns=None):
        ''' start the client
        
        Parameters
        ----------
        nruns : integer, optional
            stop after so many connect runs
        '''
        # global minimum id
        system = self.system
        pot = system.get_potential()

        self.gid = {}

        # create a fake database in memory
        db = system.create_database(db=":memory:")

        # connect to events and forward them to server
        db.on_minimum_added.connect(self._minimum_added)
        db.on_ts_added.connect(self._ts_added)

    
        while True:
            print "Obtain a new job"
            id1, coords1, id2, coords2 = self.connect_manager.get_connect_job(self.strategy)
            
            print "processing connect run between minima with global id", id1, id2
            
            print "add minima to local database"
            min1 = db.addMinimum(pot.getEnergy(coords1), coords1)
            min2 = db.addMinimum(pot.getEnergy(coords2), coords2)
            print "assigned local ids", min1._id, min2._id
            print "run double ended connect"
            connect = system.get_double_ended_connect(min1, min2, db, fresh_connect=True)
            connect.connect()
            if nruns is not None:
                nruns -= 1
                if nruns == 0: break

        print "finished sucessfully!"
        print "Data collected during run:"
        print db.session.query(Minimum).count(), "minima"
        print db.session.query(TransitionState).count(), "transition states"

    # forwarding new minima to server
    def _minimum_added(self, minimum):
        minid = self.connect_manager.add_minimum(minimum.energy, minimum.coords)
        # store the id of minimum on serverside for quick access later on
        self.gid[minimum] = minid
#        print "added minimum %d as %d"%(minimum._id, minid)
        
    # forwarding new ts to server 
    def _ts_added(self, ts):
        id1 = self.gid[ts.minimum1]
        id2 = self.gid[ts.minimum2]
        
        self.connect_manager.add_ts(id1, id2, ts.energy, ts.coords, eigenval=ts.eigenval, eigenvec=ts.eigenvec)