"""This module implements storage of pygmin simulation data in a relational database"""
import sqlalchemy
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import threading
import numpy as np
from sqlalchemy import Column, Integer, Float, PickleType
from sqlalchemy import ForeignKey
from sqlalchemy.orm import relationship, backref, deferred
from sqlalchemy.ext.declarative import declarative_base

verbose=False

Base = declarative_base()

class Minimum(Base):
    '''Base class to store minima'''
    __tablename__ = 'tbl_minima'

    _id = Column(Integer, primary_key=True)
    energy = Column(Float) 
    '''energy of the minimum'''
    # deferred means the object is loaded on demand, that saves some time / memory for huge graphs
    coords = deferred(Column(PickleType))
    '''coordinatse of the minimum'''
    
    def right_neighbors(self):
        return [x.higher_node for x in self.left_edges]

    def left_neighbors(self):
        return [x.lower_node for x in self.right_edges]

    def __init__(self, energy, coords):
        self.energy = energy
        self.coords = coords
    
#    transition_states = relationship("transition_states", order_by="transition_states.id", backref="minima")
    
class TransitionState(Base):
    '''Base class to store transition states'''
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
        self.minimum1 = min1
        self.minimum2 = min2
        self.eigenvec = eigenvec
        self.eigenval = eigenval


class Storage(object):
    '''Database storage class
    
    The Storage class handles the connection to the database.

    :Example:

    >>> from pygmin.storage import database
    >>> db = database.Storage(db="test.db")
    >>> for energy in np.random.random(10):
    >>>     a.addMinimum(energy, np.random.random(10)
    >>>
    >>> for minimum in database.minima():
    >>>     print minimum.energy
    
    '''
    engine = None
    Session = None
    accuracy = 1e-3
    onMinimumRemoved=[]
    onMinimumAdded=[]
    compareMinima=None
    
    def __init__(self, accuracy=1e-3, db=":memory:", connect_string='sqlite:///%s',\
                 onMinimumAdded=None, onMinimumRemoved=None, compareMinima=None):
        """Constructor for database storage class

        :param accuracy: energy tolerance to cound minima as equal
        :type accuracy: float
        :param db: database to connect to
        :type db: string
        :param connect_string: connection string, default is sqlite database
        :type connect_string: string
        :type message: string
        """

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
        
        - **parameters**, **types**, **return** and **return types**::
        
        :param E: energy
        :type E: float
        :param coords: coordinates of the minimum
        :type coords: np.array
        :param commit: commit changes to database
        :type commit: boolean
        :returns: minimum which was added
        :rtype: Minimum
        """
        self.lock.acquire()
        candidates = self.session.query(Minimum).\
            filter(Minimum.energy > E-self.accuracy).\
            filter(Minimum.energy<E+self.accuracy)
        
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
        return new
    
    def addTransitionState(self, E, coords, min1, min2, commit=True, eigenval=None, eigenvec=None):
        candidates = self.session.query(TransitionState).\
            filter(TransitionState.minimum1==min1).\
            filter(TransitionState.minimum1==min2).\
            filter(TransitionState.energy > E-self.accuracy).\
            filter(TransitionState.energy<E+self.accuracy)
        
        for m in candidates:
            #if(self.compareMinima):
            #    if(self.compareMinima(new, m) == False):
            #        continue
            #self.lock.release() 
            return m
        
        new = TransitionState(E, coords, min1, min2, eigenval=eigenval, eigenvec=eigenvec)
            
        self.session.add(new)
        if(commit):
            self.session.commit()
        return new
    
    def minima(self):
        '''iterate over all minima in database
        
        :returns: query for all minima in database ordered in ascendind energy
        '''
        return self.session.query(Minimum).order_by(Minimum.energy).all()
    
    def transition_states(self):
        '''iterate over all minima in database
        
        :returns: query for all minima in database ordered in ascendind energy
        '''
        return self.session.query(TransitionState).all()
    
    def minimum_adder(self):
        '''wrapper class to add minima
        
        Since pickle cannot handle pointer to member functions, this class wraps the call to
        add minimum.
        
        :returns: add minimum handler
        '''
        class minimum_adder:
            def __init__(self, db):
                self.db = db
            def __call__(self, E, coords):
                self.db.addMinimum(E, coords)
        return minimum_adder(self)

if __name__ == "__main__":    
    db = Storage()#db="test.db")
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