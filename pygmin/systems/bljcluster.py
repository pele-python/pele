from pygmin.systems import AtomicCluster
from pygmin.potentials.ljpshiftfast import LJpshift

__all__ = ["BLJCluster"]

class BLJCluster(AtomicCluster):
    """
    define the System class for a lennard jones cluster
    """
    def __init__(self, natoms, ntypeA="default", **potential_kwargs):
        self.natoms = natoms
        if ntypeA == "default":
            self.ntypeA = int(self.natoms * 0.8)
        else:
            self.ntypeA = ntypeA
        self.potential_kwargs = potential_kwargs
    
    def get_potential(self):
        return LJpshift(self.natoms, self.ntypeA, **self.potential_kwargs)
    
    def get_permlist(self):
        return [range(self.ntypeA), range(self.ntypeA, self.natoms)]
    
    def create_database(self, *args, **kwargs):
        if "accuracy" not in kwargs:
            energy_accuracy = 1e-3
            kwargs["accuracy"] = energy_accuracy
        return super(BLJCluster, self).create_database(*args, **kwargs)



def run():
    #create the system object
    sys = BLJCluster(15)
    
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