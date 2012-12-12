from pygmin.systems import AtomicCluster
from pygmin.potentials import LJ

__all__ = ["LJCluster"]

class LJCluster(AtomicCluster):
    """
    define the System class for a lennard jones cluster
    """
    def __init__(self, natoms):
        self.natoms = natoms
    
    def get_permlist(self):
        return [range(self.natoms)]
    
    def get_potential(self):
        return LJ()
    
    def create_database(self, *args, **kwargs):
        if "accuracy" not in kwargs:
            energy_accuracy = 1e-3
            kwargs["accuracy"] = energy_accuracy
        return super(LJCluster, self).create_database(*args, **kwargs)



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