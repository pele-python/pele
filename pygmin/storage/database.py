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
    __tablename__ = 'tbl_minima'

    _id = Column(Integer, primary_key=True)
    energy = Column(Float)
    # deferred means the object is loaded on demand, that saves some time / memory for huge graphs
    coords = deferred(Column(PickleType)) 
    
    def right_neighbors(self):
        return [x.higher_node for x in self.left_edges]

    def left_neighbors(self):
        return [x.lower_node for x in self.right_edges]

    def __init__(self, energy, coords):
        self.energy = energy
        self.coords = coords
    
#    transition_states = relationship("transition_states", order_by="transition_states.id", backref="minima")
    
class TransitionState(Base):
    __tablename__ = "tbl_transition_states"
    _id = Column(Integer, primary_key=True)
    energy = Column(Float)
    coords = Column(PickleType)
    
    _minimum1_id = Column(Integer, ForeignKey('tbl_minima._id'))
    minimum1 = relationship("Minimum",
                            primaryjoin="Minimum._id==TransitionState._minimum1_id",
                            backref='left_edges')

    
    _minimum2_id = Column(Integer, ForeignKey('tbl_minima._id'))
    minimum2 = relationship("Minimum",
                            primaryjoin="Minimum._id==TransitionState._minimum2_id",
                            backref='right_edges')
    
    def __init__(self, energy, min1, min2):
        self.energy = energy
        self.minimum1 = min1
        self.minimum2 = min2


class Storage(object):
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
    
    def minima(self):
        return self.session.query(Minimum).order_by(Minimum.energy).all()
    
    def minimum_adder(self):
        class minimum_adder:
            def __init__(self, db):
                self.db = db
            def __call__(self, E, coords):
                db.addMinimum(E, coords)
        return minimum_adder(self)

if __name__ == "__main__":    
    db = Storage()#db="test.db")
    m1 = db.addMinimum(1., np.random.random(10))
    db.addMinimum(3., np.random.random(10))
    db.addMinimum(2., np.random.random(10))
    db.addMinimum(1.001, np.random.random(10))
    db.minimum_adder()(7, np.random.random(10))
    
    for m in db.minima():
        print m._id, m.energy
    #for i in db.minima():
    #    print i