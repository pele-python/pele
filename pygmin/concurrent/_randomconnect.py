from random import choice
from collections import deque

import Pyro4
import sqlalchemy
import networkx as nx

from pygmin.storage import Minimum, TransitionState
from pygmin.landscape import Graph


__all__ = ["RandomConnectServer", "RandomConnectWorker"]

# we need to run pyros in multiplex mode, otherwise we run into problems with 
# SQLAlchemy. This is due to the fact that a session can only be used from one
# thread. We really should fix this issue and allow for multiple sessions!
Pyro4.config.SERVERTYPE = "multiplex"



class ConnectManagerCombine(object):
    """a class to organize which minima to try to connect
    
    Parameters
    ----------
    database : Database object
    list_len : int
        the class will create a list of minima pairs of length
        list_len.  When this list is empty the list will be rebuilt.
        Essentially this parameter indicates how often to rebuild the list
        of minima pairs to connect
    clust_min : int
        Clusters of minima below this size will be ignored.
    """
    def __init__(self, database, list_len=20, clust_min=4):
        self.db = database
        self.list_len = list_len
        self.clust_min = clust_min
        
        self.minpairs = deque()
    
    def _generate_list_combine(self):
        """make a list of minima pairs to try to connect"""
        print "analyzing the database to find minima to connect"
        self.minpairs = deque()
        
        graph = Graph(self.db).graph
        cclist = nx.connected_components(graph)

        # remove clusters with fewer than clust_min
        cclist = [cc for cc in cclist if len(cc) >= self.clust_min]
        
        if len(cclist) == 0:
            print "all minima are connected"
            return self.minpairs
        
        # get the group that all other groups will be connected to
        group1 = cclist[0]
        min1 = sorted(group1, key=lambda m: m.energy)[0]
        if True:
            # make sure that the global minimum is in group1
            global_min = self.db.minima()[0]
            if not global_min in group1:
                print "warning, the global minimum is not the in the largest cluster.  Will try to connect them"
                self.minpairs.append((min1, global_min))
                

        # remove group1 from cclist
        cclist.remove(group1)

        # get a minima from each of the other groups
        for group2 in cclist:
            if len(self.minpairs) > self.list_len:
                break

            print "adding groups of size", len(group1), "and", len(group2), "to the connect list"
            
            # sort the groups by energy
            group2.sort(key = lambda m:m.energy)
    
            # select the lowest energy minima in the groups
            # (this can probably be done in a more intelligent way)
            min2 = group2[0]

            self.minpairs.append((min1, min2))
        
        return self.minpairs
        

    def get_connect_job(self):
        if len(self.minpairs) == 0:
            self._generate_list_combine()
        if len(self.minpairs) == 0:
            return None, None
        
        min1, min2 = self.minpairs.popleft()
        return min1, min2

class ConnectManagerRandom(object):
    def __init__(self, database, Emax=None):
        self.database = database
        self.Emax = Emax
    
    def get_connect_job(self):
        """select two minima randomly"""
        query =  self.database.session.query(Minimum)
        if self.Emax is not None:
            query.filter(Minimum.energy < self.Emax)
            
        min1 = query.order_by(sqlalchemy.func.random()).first()
        min2 = query.order_by(sqlalchemy.func.random()).first()
        
        print "worker requested new job, sending minima", min1._id, min2._id
        
        return min1, min2


class ConnectManager(object):
    """class to manage which minima to try to connect
    """
    def __init__(self, database, strategy="random", list_len=20, clust_min=4, Emax=None):
        self.database = database
        self.default_strategy = strategy
        
        self.manager_random = ConnectManagerRandom(self.database, Emax)
        self.manager_combine = ConnectManagerCombine(self.database, list_len=list_len, clust_min=4)

    def get_connect_job(self, strategy=None):
        if strategy is None:
            strategy = self.default_strategy
        
        possible_strategies = ["random", "combine"]
        backup_strategy = "random"
        if strategy not in possible_strategies:
            raise Exception("strategy must be from %s" % (str(possible_strategies)))
        if strategy == "combine":
            min1, min2 = self.manager_combine.get_connect_job()
            if min1 is None or min2 is None:
                strategy = backup_strategy
            else:
                print "returning a connect job to combine two disconnected clusters"
        if strategy == "random":
            min1, min2 = self.manager_random.get_connect_job()
            print "returning a random connect job"
        
        return min1, min2

        
        

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
    
    def add_ts(self, id1, id2, E, coords):
        ''' called by worker if a new transition state is found '''
        
        print "a client found a transition state", E
        min1 = self.db.session.query(Minimum).filter(Minimum._id == id1).one()
        min2 = self.db.session.query(Minimum).filter(Minimum._id == id2).one()
        
        return self.db.addTransitionState(E, coords, min1, min2)._id

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
        print "added minimum %d as %d"%(minimum._id, minid)
        
    # forwarding new ts to server 
    def _ts_added(self, ts):
        id1 = self.gid[ts.minimum1]
        id2 = self.gid[ts.minimum2]
        
        self.connect_manager.add_ts(id1, id2, ts.energy, ts.coords)