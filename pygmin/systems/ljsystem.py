import numpy as np

from pygmin.systems import BaseSystem
from pygmin.potentials import LJ
from pygmin.takestep import RandomDisplacement, AdaptiveStepsizeTemperature
from pygmin.transition_states import orthogopt
from pygmin.mindist import minPermDistStochastic, MinDistWrapper, ExactMatchCluster


class LJCluster(BaseSystem):
    """
    define the System class for a lennard jones cluster
    """
    def __init__(self, natoms):
        self.natoms = natoms
    
    def get_potential(self):
        return LJ(self.natoms)
    
    def get_random_configuration(self):
        coords = np.random.uniform(-1, 1, 3*self.natoms) * 1.5 * float(self.natoms)**(1./3)
        return coords
    
    def get_compare_exact(self, **kwargs):
        match = ExactMatchCluster(permlist=[range(self.natoms)], **kwargs)
        return lambda m1, m2: match(m1.coords, m2.coords)
    
    def get_mindist(self, **kwargs):
        return MinDistWrapper(minPermDistStochastic, permlist=[range(self.natoms)], **kwargs)
    
    def create_database(self, *args, **kwargs):
        if "accuracy" not in kwargs:
            energy_accuracy = 1e-3
            kwargs["accuracy"] = energy_accuracy
        return super(LJCluster, self).create_database(*args, **kwargs)
        
    
    def get_takestep(self, stepsize=0.6, **kwargs):
        takeStep = RandomDisplacement(stepsize=stepsize)
        tsAdaptive = AdaptiveStepsizeTemperature(takeStep, **kwargs)
        return tsAdaptive
    
    def get_orthogonalize_to_zero_eigenvectors(self):
        return orthogopt



def run():
    #create the system object
    sys = LJCluster(15)
    
    #create a database
    db = sys.create_database()
    
    #do a short basinhopping run
    bh = sys.get_basinhopping(database=db, outstream=None)
    while len(db.minima()) < 2:
        bh.run(100)
    
    #try to connect the lowest two minima
    min1, min2 = db.minima()[:2]
    connect = sys.get_double_ended_connect(min1, min2, db)
    connect.connect()
    

if __name__ == "__main__":
    run()