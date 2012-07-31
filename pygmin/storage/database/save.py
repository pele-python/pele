import sqlalchemy
from sqlalchemy import create_engine
#from pygmin.storage.database import *
import numpy as np
from dbtypes import Minimum
import globals
from sqlalchemy.orm import sessionmaker
import threading

class SaveDB(object):
    engine = None
    Session = None
    accuracy = 1e-3
    onMinimumRemoved=[]
    onMinimumAdded=[]
    compareMinima=None
    
    def __init__(self, accuracy=1e-3, db=":memory:", connect_string='sqlite:///%s',\
                 onMinimumAdded=None, onMinimumRemoved=None, compareMinima=None):
        self.engine = create_engine(connect_string%(db), echo=True)
        globals.Base.metadata.create_all(self.engine)
        Session = sessionmaker(bind=self.engine)
        self.session = Session()
        self.accuracy=accuracy
        self.onMinimumAdded=onMinimumAdded
        self.onMinimumRemoved=onMinimumRemoved
        self.compareMinima=compareMinima
        self.lock = threading.Lock()
        
        
    def addMinimum(self, E, coords, commit=True):
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
    db = SaveDB(db="test.db")
    m1 = db.addMinimum(1., np.random.random(10))
    db.addMinimum(3., np.random.random(10))
    db.addMinimum(2., np.random.random(10))
    db.addMinimum(1.001, np.random.random(10))
    db.minimum_adder()(7, np.random.random(10))
    
    for m in db.minima():
        print m._id, m.energy
    #for i in db.minima():
    #    print i