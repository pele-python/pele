from pele.mindist import StandardClusterAlignment

__all__= ["PointGroupOrderCluster"]

class PointGroupOrderCluster(object):
    ''' Determines the point group order of a cluster
    
        Uses exact_match and standard_alignment to determine the point group 
        order of a cluster.
        
        Parameters
        ----------
        
        exact_match : 
            class to perform exact match checks and get transform policies
        
    '''
    def __init__(self, exact_match):
        self.exact_match = exact_match
        
        if not hasattr(exact_match, "transform") or not hasattr(exact_match, "measure"):
            raise RuntimeError("exact_match does not have a transform or measure policy")
        
    def __call__(self, coords):
        ''' calculate the point group order of coordinates '''
        com1 = self.exact_match.measure.get_com(coords)
        x1 = coords.copy()
        self.exact_match.transform.translate(x1, -com1)
        
        inversion_multiplier = 1
        if self.exact_match.transform.can_invert():
            x2 = x1.copy()
            self.exact_match.transform.invert(x2)   
            if self.exact_match(x1, x2, check_inversion=False):
                inversion_multiplier = 2
                
        pgorder = 0
        for rot, invert in self.exact_match.standard_alignments(x1, x1):
            if invert: continue
            if self.exact_match.check_match(x1, x1, rot, invert, reoptimize_rotation=False):
                pgorder += 1
#                print "pgorder", pgorder
        
        return inversion_multiplier * pgorder

def testlj75():
    import numpy as np
    coords = np.genfromtxt("tests/coords.lj75.gmin")
    from pele.mindist import ExactMatchAtomicCluster
    
    permlist = [range(75)]
    match = ExactMatchAtomicCluster(permlist=permlist, can_invert=True)
    calculator = PointGroupOrderCluster(match)
    pgorder = calculator(coords)
    print pgorder

def testlj6():
    from pele.systems import LJCluster
    from pele.thermodynamics import get_thermodynamic_information
    system = LJCluster(6)
    db = system.create_database()
    bh = system.get_basinhopping(db)
    bh.setPrinting(ostream=None)
    bh.run(10)
    
    get_thermodynamic_information(system, db)

    for m in db.minima():
        print m.energy, m.pgorder
        
    

if __name__ == "__main__":
    testlj6()
