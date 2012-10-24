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
from sqlalchemy.ext.declarative import declarative_base

__all__ = ["Minimum", "TransitionState", "Database"]

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
    
    Notes
    -----
    
    To avoid any double entries of minima and be able to compare them,
    only use Database.addMinimum to create a minimum object.
    
    '''
    __tablename__ = 'tbl_minima'

    _id = Column(Integer, primary_key=True)
    energy = Column(Float) 
    # deferred means the object is loaded on demand, that saves some time / memory for huge graphs
    coords = deferred(Column(PickleType))
    '''coordinates'''
    
    def __init__(self, energy, coords):
        self.energy = energy
        self.coords = np.copy(coords)
 
    def right_neighbors(self):
        return [x.higher_node for x in self.left_edges]

    def left_neighbors(self):
        return [x.lower_node for x in self.right_edges]
   
#    transition_states = relationship("transition_states", order_by="transition_states.id", backref="minima")
    
class TransitionState(Base):
    '''Transition state object 
       
    The TransitionState class represents a saddle point in the databae.     
    
    Parameters
    ----------
    energy : float
    coords : numpy array
    min1: Minimum object
        first minimum
    min2: Minimum object
        first minimum
    eigenval: float, optional
        lowest (single negative) eigenvalue of the sadle point
    eigenvec: numpy array, optional
        eigenvector which correpsonds to the negative eigenvalue 
    
    Attributes
    ----------
    energy : float
    coords : numpy array
    minimum1: Minimum object
    minimum2: Minium object
    eigenval: float
    eigenvec: numpy array
    
    Notes
    -----
    To avoid any double entries and be able to compare them, only use 
    Database.addTransitionState to create a TransitionStateobject.
    
    '''
    __tablename__ = "tbl_transition_states"
    _id = Column(Integer, primary_key=True)
    
    energy = Column(Float)
    '''energy of transition state'''
    
    coords = Column(PickleType)
    '''coordinates of transition state'''
    
    _minimum1_id = Column(Integer, ForeignKey('tbl_minima._id'))
    minimum1 = relationship("Minimum",
                            primaryjoin="Minimum._id==TransitionState._minimum1_id",
                            backref='left_edges')
    '''first minimum which connects to transition state'''
    
    _minimum2_id = Column(Integer, ForeignKey('tbl_minima._id'))
    minimum2 = relationship("Minimum",
                            primaryjoin="Minimum._id==TransitionState._minimum2_id",
                            backref='right_edges')
    '''second minimum which connects to transition state'''
    
    eigenval = Column(Float)
    '''coordinates of transition state'''

    eigenvec = Column(PickleType)
    '''coordinates of transition state'''
    
    
    def __init__(self, energy, coords, min1, min2, eigenval=None, eigenvec=None):
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
        self.dist = dist
        self.minimum1 = min1
        self.minimum2 = min2

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
    accuracy : float, optional
        energy tolerance to cound minima as equal
    db : string, optional
        database to connect to
    connect_string : string, optional
        connection string, default is sqlite database
    

    Attributes
    ----------
    engine : sqlalchemy database engine
    session : sqlalchemy session
    
    accuracy : float
    onMinimumRemoved : list of function
    onMinimumAdded : list of funtion
    compareMinima
    
    Example
    -------
    
    >>> from pygmin.storage import database
    >>> db = database.Database(db="test.db")
    >>> for energy in np.random.random(10):
    >>>     a.addMinimum(energy, np.random.random(10)
    >>>
    >>> for minimum in database.minima():
    >>>     print minimum.energy
    
    '''
    engine = None
    Session = None
    session = None
    accuracy = 1e-3
    onMinimumRemoved=[]
    onMinimumAdded=[]
    compareMinima=None
    
    def __init__(self, db=":memory:", accuracy=1e-3, connect_string='sqlite:///%s',\
                 onMinimumAdded=None, onMinimumRemoved=None, compareMinima=None):
        self.engine = create_engine(connect_string%(db), echo=verbose)
        Base.metadata.create_all(self.engine)
        Session = sessionmaker(bind=self.engine)
        self.session = Session()
        self.accuracy=accuracy
        self.onMinimumAdded=onMinimumAdded
        self.onMinimumRemoved=onMinimumRemoved
        self.compareMinima=compareMinima
        self.lock = threading.Lock()
        
        
    def addMinimum(self, E, coords, commit=True):
        """add a new minimum to database
        
        Parameters
        ----------
        energy : float
        coords : numpy.array
            coordinates of the minimum
        commit : boolean, optional
            commit changes to database, default is True
        
        Returns
        -------
        minimum : Minimum
            minimum which was added
            
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
        if(self.onMinimumAdded):
            self.onMinimumAdded(new)
        return new
    
    def addTransitionState(self, E, coords, min1, min2, commit=True, eigenval=None, eigenvec=None):
        m1, m2 = min1, min2    
        candidates = self.session.query(TransitionState).\
            filter(or_(
                       and_(TransitionState.minimum1==m1, 
                            TransitionState.minimum2==m2),
                       and_(TransitionState.minimum1==m2, 
                            TransitionState.minimum2==m1),
                       )).\
            filter(TransitionState.energy > E-self.accuracy).\
            filter(TransitionState.energy < E+self.accuracy)
        
        for m in candidates:
            #if(self.compareMinima):
            #    if(self.compareMinima(new, m) == False):
            #        continue
            #self.lock.release()
            return m

        if(m2.energy < m1.energy):
            m1,m2 = m2,m1
        new = TransitionState(E, coords, m1, m2, eigenval=eigenval, eigenvec=eigenvec)
            
        self.session.add(new)
        if(commit):
            self.session.commit()
        return new

    def getTransitionState(self, min1, min2):
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
    
    def setDistance(self, dist, min1, min2, commit=True):
        m1, m2 = min1, min2
        candidates = self.session.query(Distance).\
            filter(or_(
                       and_(Distance.minimum1==m1, 
                            Distance.minimum2==m2),
                       and_(Distance.minimum1==m2, 
                            Distance.minimum2==m1),
                       ))
        
        for m in candidates:
            #if(self.compareMinima):
            #    if(self.compareMinima(new, m) == False):
            #        continue
            #self.lock.release()
            return m
        
        new = Distance(dist, m1, m2)
            
        self.session.add(new)
        if(commit):
            self.session.commit()
        return new
    def getDistance(self, min1, min2, commit=True):
        m1, m2 = min1, min2
        candidates = self.session.query(Distance).\
            filter(or_(
                       and_(Distance.minimum1==m1, 
                            Distance.minimum2==m2),
                       and_(Distance.minimum1==m2, 
                            Distance.minimum2==m1),
                       ))

        for m in candidates:
            #if(self.compareMinima):
            #    if(self.compareMinima(new, m) == False):
            #        continue
            #self.lock.release()
            return m.dist
        return None
        
    def distances(self):
        '''iterate over all distances in database
        
        Returns
        -------
        distance : list of Distance objects
        '''
        return self.session.query(Distance).all()

        
    
    def minima(self):
        '''iterate over all minima in database
        
        Returns
        -------
        minima : list of Minimum objects
            query for all minima in database ordered in ascending energy
        '''
        return self.session.query(Minimum).order_by(Minimum.energy).all()
    
    def transition_states(self):
        '''iterate over all transition states in database
        
        Returns
        -------
        ts: list of TransitionState objects
            query for all transition states in database ordered in ascending energy
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

if __name__ == "__main__":    
    db = Database()#db="test.db")
    m1 = db.addMinimum(1., np.random.random(10))
    m2 = db.addMinimum(3., np.random.random(10))
    db.addMinimum(2., np.random.random(10))
    db.addMinimum(1.001, np.random.random(10))
    db.minimum_adder()(7, np.random.random(10))
    ts = db.addTransitionState(10., None, m1, m2)
    ts.eigenval=11.
    
    
    for m in db.minima():
        print m._id, m.energy
    for i in db.transition_states():
        print i, i.energy, i.eigenval, i.eigenvec
