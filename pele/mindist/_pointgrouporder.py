from __future__ import print_function
import numpy as np

__all__= ["PointGroupOrderCluster"]

def _rotation_in_list(rot, rot_list, eps=1e-6):
    """return true if the rotation matrix is already in the list"""
    for r in rot_list:
        if np.linalg.norm(rot - r) < eps:
            return True
    return False
            

class PointGroupOrderCluster(object):
    """ Determines the point group order of a cluster

        Uses exact_match and standard_alignment to determine the point group
        order of a cluster.

        Parameters
        ----------

        exact_match :
            class to perform exact match checks and get transform policies
    """
    def __init__(self, exact_match):
        self.exact_match = exact_match
        
        if not hasattr(exact_match, "transform") or not hasattr(exact_match, "measure"):
            raise RuntimeError("exact_match does not have a transform or measure policy")
        
    def __call__(self, coords):
        """ calculate the point group order of coordinates """
        com1 = self.exact_match.measure.get_com(coords)
        x1 = coords.copy()
        self.exact_match.transform.translate(x1, -com1)
        
        inversion_multiplier = 1
        if self.exact_match.transform.can_invert():
            x2 = x1.copy()
            self.exact_match.transform.invert(x2)   
            if self.exact_match(x1, x2, check_inversion=False):
                inversion_multiplier = 2
                
        rot_list = []
        for rot, invert in self.exact_match.standard_alignments(x1, x1):
            if invert: continue
            match = self.exact_match.check_match(x1, x1, rot, invert)
            if match is not None:
                # check if the rotation matrix has already been found.
                if not _rotation_in_list(match.rotation, rot_list):
                    rot_list.append(match.rotation)

        pgorder = len(rot_list)
        return inversion_multiplier * pgorder

#
# testing only below here
#

def testlj75(): # pragma: no cover
    import numpy as np
    coords = np.genfromtxt("tests/coords.lj75.gmin")
    from pele.mindist import ExactMatchAtomicCluster
    
    permlist = [list(range(75))]
    match = ExactMatchAtomicCluster(permlist=permlist, can_invert=True)
    calculator = PointGroupOrderCluster(match)
    pgorder = calculator(coords)
    print(pgorder)

def testlj6(): # pragma: no cover
    from pele.systems import LJCluster
    from pele.thermodynamics import get_thermodynamic_information
    system = LJCluster(6)
    db = system.create_database()
    bh = system.get_basinhopping(db)
    bh.setPrinting(ostream=None)
    bh.run(10)
    
    get_thermodynamic_information(system, db)

    for m in db.minima():
        print(m.energy, m.pgorder)
        
    

if __name__ == "__main__":
    testlj6()

