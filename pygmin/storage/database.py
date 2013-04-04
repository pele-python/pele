"""Database for simulation data in a relational database
"""
import sqlalchemy
from sqlalchemy import create_engine, and_, or_
from sqlalchemy.orm import sessionmaker
import threading
import numpy as np
from sqlalchemy import Column, Integer, Float, PickleType
from sqlalchemy import ForeignKey
from sqlalchemy.orm import relationship, backref, deferred
import sqlalchemy.orm
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.sql import select, bindparam, case, insert
from sqlalchemy.schema import Index
from pygmin.utils.events import Signal
import os

__all__ = ["Minimum", "TransitionState", "Database", "Distance"]

_schema_version = 1
verbose=False

Base = declarative_base()

class Minimum(Base):
    '''    
    The Minimum class represents a minimum in the database.
    
    Parameters
    ----------
    energy : float
    coords : numpy array
        coordinates
    
    Attributes
    ----------
    energy : float
    coords : numpy array
    fvib : float
        log product of squared frequencies for free energy calculation
    pgorder : integer
        point group order
        
    Notes
    -----
    
    To avoid any double entries of minima and be able to compare them,
    only use `Database.addMinimum()` to create a minimum object.

    See Also
    --------
    Database, TransitionState, Distance
    '''
    __tablename__ = 'tbl_minima'

    _id = Column(Integer, primary_key=True)
    energy = Column(Float) 
    # deferred means the object is loaded on demand, that saves some time / memory for huge graphs
    coords = deferred(Column(PickleType))
    fvib = Column(Float)
    pgorder = Column(Integer)
    
    '''coordinates'''
    
    def __init__(self, energy, coords):
        self.energy = energy
        self.coords = np.copy(coords)
 
    def right_neighbors(self):
        return [x.higher_node for x in self.left_edges]

    def left_neighbors(self):
        return [x.lower_node for x in self.right_edges]
   
    def __eq__(self, m):
        """m can be integer or Minima object"""
        assert self._id is not None
        if isinstance(m, Minimum):
            assert m._id is not None
            return self._id == m._id
        else:
            return self._id == m
        
    def __hash__(self):
        assert self._id is not None
        return self._id
         
#    transition_states = relationship("transition_states", order_by="transition_states.id", backref="minima")
    
class TransitionState(Base):
    '''Transition state object 
       
    The TransitionState class represents a saddle point in the database.     
    
    Parameters
    ----------
    energy : float
    coords : numpy array
    min1 : Minimum object
        first minimum
    min2 : Minimum object
        first minimum
    eigenval : float, optional
        lowest (single negative) eigenvalue of the saddle point
    eigenvec : numpy array, optional
        eigenvector which corresponds to the negative eigenvalue 
    fvib : float
        log product of squared frequencies for free energy calculation
    pgorder : integer
        point group order
    
    
    
    Attributes
    ----------
    energy : float
    coords : numpy array
    minimum1 : Minimum object
    minimum2 : Minium object
    eigenval : float
    eigenvec : numpy array
    
    Notes
    -----
    To avoid any double entries and be able to compare them, only use 
    Database.addTransitionState to create a TransitionStateobject.
    
    programming note: The functions in the database require that 
    ts.minimum1._id < ts.minimum2._id.  This will be handled automatically
    by the database, but we must remember not to screw it up

    See Also
    --------
    Database, Minimum, Distance
    '''
    __tablename__ = "tbl_transition_states"
    _id = Column(Integer, primary_key=True)
    
    energy = Column(Float)
    '''energy of transition state'''
    
    coords = deferred(Column(PickleType))
    '''coordinates of transition state'''
    
    _minimum1_id = Column(Integer, ForeignKey('tbl_minima._id'))
    minimum1 = relationship("Minimum",
                            primaryjoin="Minimum._id==TransitionState._minimum1_id") #,
                            #backref='left_edges')
    '''first minimum which connects to transition state'''
    
    _minimum2_id = Column(Integer, ForeignKey('tbl_minima._id'))
    minimum2 = relationship("Minimum",
                            primaryjoin="Minimum._id==TransitionState._minimum2_id")#,
                            #backref='right_edges')
    '''second minimum which connects to transition state'''
    
    eigenval = Column(Float)
    '''coordinates of transition state'''

    eigenvec = deferred(Column(PickleType))
    '''coordinates of transition state'''

    fvib = Column(Float)
    pgorder = Column(Integer)    
    
    def __init__(self, energy, coords, min1, min2, eigenval=None, eigenvec=None):
        assert min1._id is not None
        assert min2._id is not None
        
        self.energy = energy
        self.coords = np.copy(coords)
        if(min1._id < min2._id):
            self.minimum1 = min1
            self.minimum2 = min2
        else:
            self.minimum1 = min2
            self.minimum2 = min1
            
        self.eigenvec = np.copy(eigenvec)
        self.eigenval = eigenval

class Distance(Base):
    '''object to store "mindist" distances between minima
    
    Parameters
    ----------
    dist: float
        distance between minima
    min1: Minimum object
        first minimum
    min2: Minimum object
        second minimum
        
    Attributes
    ----------
    dist: float
    minimum1: Minimum object
    minimum2: Minimum object
    
    Notes
    -----
    To avoid any double entries and be able to compare them,
    only use `Database.setDistance()` to create a Distance object.
    
    programming note: As with TransitionState, functions in the database require that 
    `dist.minimum1._id < dist.minimum2._id`.  This will be handled automatically
    by the database, but we must remember not to screw it up.
    
    See Also
    --------
    Database, Minimum, TransitionState

    '''
    __tablename__ = "tbl_distances"
    _id = Column(Integer, primary_key=True)
    
    dist = Column(Float)
    '''distance between minima'''
        
    _minimum1_id = Column(Integer, ForeignKey('tbl_minima._id'))
    minimum1 = relationship("Minimum",
                            primaryjoin="Minimum._id==Distance._minimum1_id")
    '''first minimum'''
    
    _minimum2_id = Column(Integer, ForeignKey('tbl_minima._id'))
    minimum2 = relationship("Minimum",
                            primaryjoin="Minimum._id==Distance._minimum2_id")
    '''second minimum'''
    
    
    def __init__(self, dist, min1, min2):
        assert min1._id is not None
        assert min2._id is not None
        
        if(min1._id < min2._id):
            m1, m2 = min1, min2
        else:
            m1, m2 = min2, min1
            
        self.dist = dist
        self.minimum1 = min1
        self.minimum2 = min2

Index('idx_transition_states', TransitionState.__table__.c._minimum1_id, TransitionState.__table__.c._minimum2_id)
Index('idx_distances', Distance.__table__.c._minimum1_id, Distance.__table__.c._minimum2_id, unique=True)


class Database(object):
    '''Database storage class
    
    The Database class handles the connection to the database. It has functions to create new Minima,
    TransitionState and Distance objects. The objects are persistent in the database and exist as
    soon as the Database class in connected to the database. If any value in the objects is changed,
    the changes are automatically persistent in the database (TODO: be careful, check commit transactions, ...)
    
    Database uses SQLAlchemy to connect to the database. Check the web page for available connectors. Unless
    you know better, the standard sqlite should be used. The database can be generated in memory (default) or
    written to a file if db is specified when creating the class.

    Parameters
    ----------
    db : string, optional
        filename of new or existing database to connect to.  default creates
        new database in memory.
    accuracy : float, optional
        energy tolerance to count minima as equal
    connect_string : string, optional
        connection string, default is sqlite database
    compareMinima : callable, `bool = compareMinima(min1, min2)`, optional
        called to determine if two minima are identical.  Only called
        if the energies are within `accuracy` of each other.
    createdb : boolean, optional
        create database if not exists, default is true
        
    Attributes
    ----------
    engine : sqlalchemy database engine
    session : sqlalchemy session
    
    accuracy : float
    on_minimum_removed : signal 
        called when a minimum is removed from the database 
    on_minimum_added : signal
        called when a new, unique, minimum is added to the database
    on_ts_removed : signal 
        called when a transition_state is removed from the database 
    on_ts_added : signal
        called when a new, unique, transition state is added to the database
    compareMinima
    
    Examples
    --------
    
    >>> from pygmin.storage import Database
    >>> db = Database(db="test.db")
    >>> for energy in np.random.random(10):
    >>>     a.addMinimum(energy, np.random.random(10))
    >>>
    >>> for minimum in database.minima():
    >>>     print minimum.energy
    
    See Also
    --------
    Minimum
    TransitionState
    Distance
    
    '''
    engine = None
    Session = None
    session = None
    connection = None
    accuracy = 1e-3
    compareMinima=None
        
    def __init__(self, db=":memory:", accuracy=1e-3, connect_string='sqlite:///%s',
                 compareMinima=None, createdb=True):
        global _schema_version
        if not createdb:
            if not os.path.isfile(db): 
                raise IOError("database does not exist")
            
        self.engine = create_engine(connect_string%(db), echo=verbose)
        if createdb:
            conn = self.engine.connect()
            if not self.engine.has_table("tbl_minima"):
                conn.execute("PRAGMA user_version = %d;"%_schema_version)
            Base.metadata.create_all(self.engine)
            result=conn.execute("PRAGMA user_version;")
            schema = result.fetchone()[0]
            result.close()
            conn.close()
            if _schema_version != schema:
                raise IOError("database schema outdated, current (newest) version: "
                              "%d (%d). Please use migrate_db.py in pygmin/scripts to update database"%(schema, _schema_version))
            
        Session = sessionmaker(bind=self.engine)
        self.session = Session()
        self.accuracy=accuracy
        self.on_minimum_added = Signal()
        self.on_minimum_removed = Signal()
        self.on_ts_added = Signal()
        self.on_ts_removed = Signal()
        
        self.compareMinima=compareMinima
        self.lock = threading.Lock()
        self.connection = self.engine.connect()
        
        self._initialize_queries()
        
    def _initialize_queries(self):
        #        self._sql_get_dist = select([Distance.__table__.c.dist],
        #               or_(and_(Distance.__table__.c._minimum1_id==bindparam("id1"), 
        #                        Distance.__table__.c._minimum2_id==bindparam("id2")),
        #                   and_(Distance.__table__.c._minimum1_id==bindparam("id2"), 
        #                        Distance.__table__.c._minimum2_id==bindparam("id1")),
        #               ), use_labels=False).limit(1)
        tbl =  Distance.__table__.c
        self._sql_get_dist = select([tbl.dist],and_(
                                 tbl._minimum1_id==bindparam("id1"), 
                                 tbl._minimum2_id==bindparam("id2")
                                 ), use_labels=False).limit(1)
                                 
        #self._sql_set_dist = Distance.__table__.insert().values(_minimum1_id=bindparam("id1"),_minimum2_id=bindparam("id2"), dist=bindparam("dist"))
        self._sql_set_dist = "INSERT OR REPLACE INTO tbl_distances (dist, _minimum1_id, _minimum2_id) VALUES (:dist, :id1, :id2)"
        
        #self._sql_set_dist_upd = Distance.__table__.update().where(and_(tbl._minimum1_id==bindparam("id1"),tbl._minimum2_id==bindparam("id2"))).values(dist=bindparam("dist"))
        #self._sql_set_dist_ins = Distance.__table__.insert().values(_minimum1_id=bindparam("id1"),_minimum2_id=bindparam("id2"), dist=bindparam("dist"))
        
    
    def addMinimum(self, E, coords, commit=True):
        """add a new minimum to database
        
        Parameters
        ----------
        energy : float
        coords : numpy.array
            coordinates of the minimum
        commit : bool, optional
            commit changes to database
        
        Returns
        -------
        minimum : Minimum
            minimum which was added (not necessarily a new minimum)
            
        """
        self.lock.acquire()
        candidates = self.session.query(Minimum).\
            filter(Minimum.energy > E-self.accuracy).\
            filter(Minimum.energy < E+self.accuracy)
        
        new = Minimum(E, coords)
            
        for m in candidates:
            if(self.compareMinima):
                if(self.compareMinima(new, m) == False):
                    continue
            self.lock.release() 
            return m
        self.session.add(new)
        if(commit):
            self.session.commit()
        self.lock.release()
        self.on_minimum_added(new)
        return new
    
    def getMinimum(self, id):
        """return the minimum with a given id"""
        return self.session.query(Minimum).get(id)
        
    def addTransitionState(self, energy, coords, min1, min2, commit=True, eigenval=None, eigenvec=None):
        """Add transition state object
        
        Parameters
        ----------
        energy : float
            energy of transition state
        coords : numpy array
            coordinates of transition state
        min1, min2 : Minimum
            minima on either side of the transition states
        eigenval : float
            the eigenvalue (curvature) across the transition state
        eigenvec : numpy array
            the eigenvector associated with eigenval
        commit : bool
            commit changes to sql database
        
        Returns
        -------
        ts : TransitionState
            the transition state object (not necessarily new)
        """
        m1, m2 = min1, min2
        if m1._id > m2._id:
            m1, m2 = m2, m1
        candidates = self.session.query(TransitionState).\
            filter(or_(
                       and_(TransitionState.minimum1==m1, 
                            TransitionState.minimum2==m2),
                       and_(TransitionState.minimum1==m2, 
                            TransitionState.minimum2==m1),
                       )).\
            filter(TransitionState.energy > energy-self.accuracy).\
            filter(TransitionState.energy < energy+self.accuracy)
        
        for m in candidates:
            #if(self.compareMinima):
            #    if(self.compareMinima(new, m) == False):
            #        continue
            #self.lock.release()
            return m

        #if(m2.energy < m1.energy):
        #    m1,m2 = m2,m1
        new = TransitionState(energy, coords, m1, m2, eigenval=eigenval, eigenvec=eigenvec)
            
        self.session.add(new)
        if(commit):
            self.session.commit()
        self.on_ts_added(new)
        return new

    def getTransitionState(self, min1, min2):
        """return the TransitionState between two minima
        
        Returns
        -------
        ts : None or TransitionState
        """
        m1, m2 = min1, min2
        candidates = self.session.query(TransitionState).\
            filter(or_(
                       and_(TransitionState.minimum1==m1, 
                            TransitionState.minimum2==m2),
                       and_(TransitionState.minimum1==m2, 
                            TransitionState.minimum2==m1),
                       ))

        for m in candidates:
            #if(self.compareMinima):
            #    if(self.compareMinima(new, m) == False):
            #        continue
            #self.lock.release()
            return m
        return None
    
    def setDistance(self, dist, min1, min2):
        """set the distance between two minima
        """
        id1 = max(min1._id, min2._id)
        id2 = min(min1._id, min2._id)
        
        self.connection.execute(self._sql_set_dist, [{'id1':id1, 'id2':id2, 'dist':dist}])
        
        #res = self.connection.execute(self._sql_set_dist_upd, [{'id1':min1._id, 'id2':min2._id, 'dist':dist}])
        #if(res.rowcount == 0):
        #    self.connection.execute(self._sql_set_dist_ins, [{'id1':min1._id, 'id2':min2._id, 'dist':dist}])
        
        
    def setDistanceBulk(self, values):
        """set multiple distances
        
        This is faster than calling `setDistance` multiple times
        
        Parameters
        -----------
        values : iterable of tuples of form ((min1, min2), dist)
        """
        submit = []
        for mins, dist in values:
            submit.append({'id1':min(mins[0]._id, mins[1]._id), 'id2':max(mins[0]._id, mins[1]._id), 'dist':dist})
        self.connection.execute(self._sql_set_dist, submit)
        
#    def setDistanceMultiple(self, newdistances, commit=True):
#        """set multiple distances at once
#        
#        this should be much faster than calling setDistance
#        multiple times.  Especially for databases with many 
#        distances
#        
#        Parameters
#        ----------
#        distances :
#            a dictionary of distances with the key a tuple of minima
#        """
#        #copy distances so it can be safely modified
#        newdistances = newdistances.copy()
#        
#        #update distances that are already in the database and remove
#        #them from newdistances
#        for d in self.distances():
#            m1id, m2id = d._minimum1_id, d._minimum2_id
#            
#            dnew = newdistances.get((m1id, m2id))
#            if dnew is not None:
#                d.dist = dnew
#                newdistances.pop((m1id, m2id))
#
#            #check alternate ordering as well
#            dnew = newdistances.get((m2id, m1id))
#            if dnew is not None:
#                d.dist = dnew
#                newdistances.pop((m2id, m1id))
#        
#        #all the rest of the distances in newdistances are not in the database.
#        #Add them all.
#        for d in newdistances.items():
#            self.session.add(Distance(d[1], d[0][0], d[0][1]))
#        
#        if(commit):
#            self.session.commit()        
            
    def getDistanceORM(self, min1, min2):
        """slow way of getting distance.  deprecated"""
        if(min1._id > min2._id):
            min1, min2 = min2, min1
            
        candidates = self.session.query(Distance).\
            filter(and_(Distance.minimum1==min1, 
                            Distance.minimum2==min2))
        try:
            return candidates.one().dist        
        except sqlalchemy.orm.exc.NoResultFound:
            return None
    
    def getDistance(self, min1, min2):
        """return the distance between two minima
        
        Returns
        --------
        dist : float or None
        """        
        result = self.connection.execute(self._sql_get_dist, id1=min1._id, id2=min2._id)
        dist = result.fetchone()
        result.close()
        if dist is None:
            return None
        return dist[0]
        
    def distances(self):
        '''return an iterator over all distances in database
        '''
        #return self.session.query(Distance).all()
        return self.session.query(Distance).yield_per(100)
        #return self.session.query(Distance)

        
    
    def minima(self, order_energy=True):
        '''return an iterator over all minima in database
        
        Parameters
        ----------
        order_energy : bool
            order the minima by energy
        '''
        if order_energy:
            return self.session.query(Minimum).order_by(Minimum.energy).all()
        else:
            return self.session.query(Minimum).all()
    
    def transition_states(self):
        '''return an iterator over all transition states in database
        '''
        return self.session.query(TransitionState).all()
    
    def minimum_adder(self, Ecut=None):
        '''wrapper class to add minima
        
        Since pickle cannot handle pointer to member functions, this class wraps the call to
        add minimum.
        
        Parameters
        ----------
        Ecut: float, optional
             energy cutoff, don't add minima which are higher in energy
        
        Returns
        -------
        handler: minimum_adder class
            minimum handler to add minima
            
        '''
        class minimum_adder:
            def __init__(self, db, Ecut):
                self.db = db
                self.Ecut = Ecut
            def __call__(self, E, coords):
                if(not self.Ecut is None):
                    if(E > self.Ecut):
                        return None
                self.db.addMinimum(E, coords)
        return minimum_adder(self, Ecut)
    
    def removeMinimum(self, m):
        """remove a minimum from the database
        
        Remove a minimum and any objects (TransitionState or Distance) 
        pointing to that minimum.
        """
        #delete any distance objects pointing to min2
        candidates = self.session.query(Distance).\
            filter(or_(Distance.minimum1 == m, 
                       Distance.minimum2 == m))
        candidates = list(candidates)
        for d in candidates:
            self.session.delete(d)
            
        #delete any transition states objects pointing to min2
        candidates = self.session.query(TransitionState).\
            filter(or_(TransitionState.minimum1 == m, 
                       TransitionState.minimum2 == m))
        candidates = list(candidates)
        for ts in candidates:
            self.on_ts_removed(ts)
            self.session.delete(ts)
        
        self.on_minimum_removed(m)
        #delete the minimum
        self.session.delete(m)
        self.session.commit()

    
    def mergeMinima(self, min1, min2):
        """merge two minima in the database
        
        min2 will be deleted and everything that 
        points to min2 will point to min1 instead.
        
        (actually, any Distance objects pointing to 
        min2 will just be deleted)
        """
        #find all transition states for which ts.minimum1 is min2
        candidates = self.session.query(TransitionState).\
            filter(TransitionState.minimum1 == min2) 
        for ts in candidates:
            #should we check if this will duplicate an existing transition state?
            ts.minimum1 = min1
            if ts.minimum1._id > ts.minimum2._id:
                ts.minimum1, ts.minimum2 = ts.minimum2, ts.minimum1
        
        #find all transition states for which ts.minimum2 is min2
        candidates = self.session.query(TransitionState).\
            filter(TransitionState.minimum2 == min2) 
        for ts in candidates:
            #should we check if this will duplicate an existing transition state?
            ts.minimum2 = min1
            if ts.minimum1._id > ts.minimum2._id:
                ts.minimum1, ts.minimum2 = ts.minimum2, ts.minimum1
              
        
        #delete any distance objects pointing to min2
        candidates = self.session.query(Distance).\
            filter(or_(Distance.minimum1 == min2, 
                       Distance.minimum2 == min2))

        #copy it into list format so the iterator doesn't get corrupted as we delete things            
        candidates = list(candidates)
        for d in candidates:
            self.session.delete(d)
        
        self.session.delete(min2)
        self.session.commit()
    
    def number_of_minima(self):
        """return the number of minima in the database
        
        Notes
        -----
        This is much faster than len(database.minima()), but is is not instantaneous.  
        It takes a longer time for larger databases.  The first call to number_of_minima() 
        can be much faster than subsequent calls.  
        """
        return self.session.query(Minimum).count()

    def number_of_transition_states(self):
        """return the number of transition states in the database
        
        Notes
        -----
        see notes for number_of_minima()
        
        See Also
        --------
        number_of_minima
        """
        return self.session.query(TransitionState).count()


if __name__ == "__main__":    
    db = Database()
    m1 = db.addMinimum(1., np.random.random(10))
    m2 = db.addMinimum(3., np.random.random(10))
    db.addMinimum(2., np.random.random(10))
    db.addMinimum(1.001, np.random.random(10))
    db.minimum_adder()(7, np.random.random(10))
    ts = db.addTransitionState(10., None, m1, m2)
    ts.eigenval=11.
    
    db.setDistance(21., m2, m1)
    db.setDistance(12., m1, m2)
    db.setDistance(211., m2, m1)
    db.setDistance(11., m1, m1)
    db.setDistance(22., m2, m2)
    #db.setDistance(11., m1, m2)
    print db.getDistance(m1, m2)
    print "List of all distances"
    for d in db.distances():
        print d.minimum1._id, d.minimum2._id, d.dist        

    print    
    for m in db.minima():
        print m._id, m.energy
    #for i in db.transition_states():
    #    print i, i.energy, i.eigenval, i.eigenvec
