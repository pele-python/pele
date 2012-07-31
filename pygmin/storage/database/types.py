from sqlalchemy import Column, Integer, Float, PickleType
from sqlalchemy import ForeignKey
from sqlalchemy.orm import relationship, backref
import globals

class Minimum(globals.Base):
    __tablename__ = 'minima'

    id = Column(Integer, primary_key=True)
    energy = Column(Float)
    coords = Column(PickleType)
    
    def __init__(self, energy, coords):
        self.energy = energy
        self.coords = coords
    
#    transition_states = relationship("transition_states", order_by="transition_states.id", backref="minima")
    
class TransitionState(globals.Base):
    __tablename__ = "transition_states"
    id = Column(Integer, primary_key=True)
    energy = Column(Float)
    coords = Column(PickleType)
    
    _minimum1_id = Column(Integer, ForeignKey('minima.id'))
    minimum1 = relationship("Minimum")

    #__minimum2_id = Column(Integer, ForeignKey('minima.id'))
    #minimum2 = relationship("Minimum")
    
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
    session.flush()
    ts = TransitionState(3., min1, min2)
    session.add(TransitionState(3., min1, min2))
    
    session.flush()
    print min1.id, min2.id, ts.id, ts._minimum1_id
    print min1.coords