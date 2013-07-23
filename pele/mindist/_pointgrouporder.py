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
            if self.exact_match.check_match(x1, x1, rot, invert):
                pgorder += 1
        
        return inversion_multiplier*pgorder