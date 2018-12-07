from __future__ import print_function
from . import Pyro4

from pele.landscape import ConnectManager


__all__ = ["ConnectServer", "ConnectWorker", "BasinhoppingWorker"]

# we need to run pyros in multiplex mode, otherwise we run into problems with 
# SQLAlchemy. This is due to the fact that a session can only be used from one
# thread. We really should fix this issue and allow for multiple sessions!
Pyro4.config.SERVERTYPE = "multiplex"
        

class ConnectServer(object):
    """
    Server which receives requests from, and passes connect jobs to the workers

    The server also receives minima and transition states from the workers
    and adds them to the database.

    Parameters
    ----------
    system : pele.system.BaseSystem
        system class to process
    database : pele.storage.Database
        working database
    server_name : string, optional
        Unique name for clients to connect to this server on current host
        (objid for pyros). None for random
    host : string, optional
        host to setup server. default is localhost which does not allow connections
        from remote machines
    port : integer, optional
        port to listen for connections

    See Also
    --------
    ConnectWorker
    pele.landscape.ConnectManager
    """
    
    def __init__(self, system, database, server_name=None, host=None, port=0):
        self.system = system
        self.db = database
        self.server_name = server_name
        self.host=host
        self.port=port
        
        self.connect_manager = ConnectManager(self.db)

    def set_connect_manager(self, connect_manager):
        """add a custom connect manager
        
        the connect manager decides which connect jobs should be performed
        """
        self.connect_manager = connect_manager
        
#    def set_emax(self, Emax):
#        raise Exception("set_emax is not implemented yet in the new ConnectManager scheme")
#        self.Emax = None

    def get_connect_job(self, strategy="random"):
        """ get a new connect job """
        min1, min2 = self.connect_manager.get_connect_job(strategy)
        return min1.id(), min1.coords, min2.id(), min2.coords

    def get_system(self):
        """ provide system class to worker """
        return self.system
    
    def add_minimum(self, E, coords):
        """ called by worker if a new minimum is found

        Returns
        -------
        ID : global id of minimum added.
        """
        print("a client found a minimum", E)
        m = self.db.addMinimum(E, coords)
        return m.id()
    
    def add_ts(self, id1, id2, E, coords, eigenval=None, eigenvec=None):
        """called by worker if a new transition state is found

        Parameters
        ----------
        id1, id2 : int
            The server-side (global) ID's of the minima on either side of the transition state.
            The worker is responsible for knowing the global id of the minima.  This ID is returned
            to the worker when a minimum is added
        E : float
            energy of transition state
        coords : array
            coordinates of transition state

        Returns
        -------
        ID : global id of transition state added
        """
        print("a client found a transition state", E)
        min1 = self.db.getMinimum(id1)
        min2 = self.db.getMinimum(id2)
        
        ts = self.db.addTransitionState(E, coords, min1, min2, eigenval=eigenval, eigenvec=eigenvec)
        return ts.id()

    def run(self):
        """ start the server and listen for incoming connections """
        print("Starting Pyros daemon")
        daemon=Pyro4.Daemon(host=self.host, port=self.port)
        # make the connect_server available to Pyros children
        uri=daemon.register(self, objectId=self.server_name)
        print("The connect server can be accessed by the following uri: ", uri)
        
        print("Ready to accept connections")
        daemon.requestLoop() 
        
class ConnectWorker(object):
    """
    worker class to execute connect runs.

    The worker will return all minima and transiton states found to the server

    Parameters
    ----------
    uri : string
        uri for job server
    system : BaseSystem, optional
        if no system class is specified, the worker obtains the system
        class from the ConnectServer by get_system. This only works for pickleable
        systems classes. If this is not the case, the system class can be
        created on the client side and passed as a parameter.
    strategy : str
        strategy to use when choosing which minima to connect

    See Also
    --------
    ConnectServer
    pele.landscape.ConnectManager
    """
    
    def __init__(self, uri, system=None, strategy="random"):
        print("connecting to",uri)
        self.connect_server = Pyro4.Proxy(uri)
        if system is None:
            system = self.connect_server.get_system()
        self.system = system
        
        self.strategy = strategy
        
    def run(self, nruns=None):
        """ start the client

        Parameters
        ----------
        nruns : integer, optional
            stop after so many connect runs
        """
        # global minimum id
        system = self.system
        pot = system.get_potential()

        # stores the global id's of the minima found
        self.gid = dict()

        # create a local database in memory
        db = system.create_database(db=":memory:")

        # connect to events and forward them to server
        db.on_minimum_added.connect(self._minimum_added)
        db.on_ts_added.connect(self._ts_added)
    
        while True:
            print("Obtain a new job")
            id1, coords1, id2, coords2 = self.connect_server.get_connect_job(self.strategy)
            
            print("processing connect run between minima with global id", id1, id2)
            
            # add minima to local database
            min1 = db.addMinimum(pot.getEnergy(coords1), coords1)
            min2 = db.addMinimum(pot.getEnergy(coords2), coords2)
            # assigned local ids", min1.id(), min2.id()

            # run double ended connect
            connect = system.get_double_ended_connect(min1, min2, db, fresh_connect=True)
            connect.connect()
            if nruns is not None:
                nruns -= 1
                if nruns == 0: break

        print("finished successfully!")
        print("Data collected during run:")
        print(db.number_of_minima(), "minima")
        print(db.number_of_transition_states(), "transition states")

    def _minimum_added(self, minimum):
        """forward new minimum to server"""
        minid = self.connect_server.add_minimum(minimum.energy, minimum.coords)
        # store the id of minimum on server-side for quick access later on
        self.gid[minimum] = minid
        
    def _ts_added(self, ts):
        """forward new transition state to server""" 
        id1 = self.gid[ts.minimum1]
        id2 = self.gid[ts.minimum2]
        
        self.connect_server.add_ts(id1, id2, ts.energy, ts.coords, eigenval=ts.eigenval, eigenvec=ts.eigenvec)

class BasinhoppingWorker(object):
    """
    worker class to execute basinhopping runs in parallel

    The worker will return all new minima to the server.

    The basinhopping run will be set up using system.get_basinhopping().

    Parameters
    ----------
    uri : string
        uri for job server
    system : BaseSystem, optional
        if no system class is specified, the worker obtains the system
        class from the ConnectServer by get_system. This only works for pickleable
        systems classes. If this is not the case, the system class can be
        created on the client side and passed as a parameter.

    See Also
    --------
    ConnectServer
    pele.Basinhopping
    """
    
    def __init__(self,uri, system=None, **basinhopping_kwargs):
        print("connecting to",uri)
        self.connect_server = Pyro4.Proxy(uri)
        if system is None:
            system = self.connect_server.get_system()
        self.system = system
        self.basinhopping_kwargs = basinhopping_kwargs
    
    def run(self, nsteps=10000):
        """ start the client

        Parameters
        ----------
        nsteps : integer
            number of basinhopping iterations
        """
        # create a local database in memory
        db = self.system.create_database(db=":memory:", **self.basinhopping_kwargs)

        # connect to events and forward them to server
        db.on_minimum_added.connect(self._minimum_added)

        bh = self.system.get_basinhopping(database=db)
        bh.run(nsteps)

        print("finished successfully!")
        print("minima found:", db.number_of_minima())

    def _minimum_added(self, minimum):
        """forward new minimum to server"""
        self.connect_server.add_minimum(minimum.energy, minimum.coords)
   
