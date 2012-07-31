from sqlalchemy import Column, Integer, Float, PickleType
from sqlalchemy import ForeignKey
from sqlalchemy.orm import relationship, backref, deferred
import globals

class Minimum(globals.Base):
    __tablename__ = 'tbl_minima'

    _id = Column(Integer, primary_key=True)
    energy = Column(Float)
    # deferred means the object is loaded on demand, that saves some time / memory for huge graphs
    coords = deferred(Column(PickleType)) 
    
    def __init__(self, energy, coords):
        self.energy = energy
        self.coords = coords
    
#    transition_states = relationship("transition_states", order_by="transition_states.id", backref="minima")
    
class TransitionState(globals.Base):
    __tablename__ = "tbl_transition_states"
    _id = Column(Integer, primary_key=True)
    energy = Column(Float)
    coords = Column(PickleType)
    
    _minimum1_id = Column(Integer, ForeignKey('tbl_minima._id'))
    minimum1 = relationship("Minimum", primaryjoin="Minimum._id==TransitionState._minimum1_id")

    
    _minimum2_id = Column(Integer, ForeignKey('tbl_minima._id'))
    minimum2 = relationship("Minimum", primaryjoin="Minimum._id==TransitionState._minimum2_id")
    
    def __init__(self, energy, min1, min2):
        self.energy = energy
        self.minimum1 = min1
        self.minimum2 = min2

if __name__ == "__main__":
    import numpy as np
    from pygmin.storage import database as db

    import sqlalchemy
    from sqlalchemy import create_engine
    
    engine = create_engine('sqlite:///:memory:', echo=True)
    
    db.Base.metadata.create_all(engine)
    from sqlalchemy.orm import sessionmaker
    Session = sessionmaker(bind=engine)
    session = Session()
    
    min1 = Minimum(1., np.random.random(10))
    min2 = Minimum(2., np.random.random(10))
    
    session.add_all([min1, min2])
    #session.flush()
    #ts = TransitionState(3., min1, min2)
    #session.add(ts)
    
    #session.flush()
    #print min1._id, min2._id, ts._id, ts._minimum1_id
    #print min1.coords