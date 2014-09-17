from potential import potential

__all__ = ["ColdFusionCheck", "addColdFusionCheck"]

class ColdFusionCheck(potential):
    """fixes for potentials that have cold fusion
    
    some potentials experience cold fusion, is our term for when the energy
    becomes un-physically small.  This is usually when two atoms sit on top of each other. 
    """
    def __init__(self, potential, coldfusionlimit=-1000, coldfusionenergy=1000000):
        self.potential = potential
        self.coldfusionlimit = coldfusionlimit
        self.coldfusionenergy = coldfusionenergy
        
    def getEnergy(self, coords):
        E = self.potential.getEnergy(coords)
        if E < self.coldfusionlimit:
            return self.coldfusionenergy
        return E
    
    def getEnergyGradient(self, coords):
        E, grad = self.potential.getEnergyGradient(coords)
        if E < self.coldfusionlimit:
            E = self.coldfusionenergy
            grad[:] = 0.
        return E, grad
    
    def __call__(self, E, conf, **kwargs):
        """test whether this energy satisfies the cold fusion criterion"""
        if E < self.coldfusionlimit:
            print "Cold fusion detected, energy was " + str(E)
            #kwargs['driver'].trial_energy = self.coldfusionenergy
            return False
        return True
        
def addColdFusionCheck(basinhopping, coldfusionlimit=-1000, coldfusionenergy=1000000):
    pot = basinhopping.potential
    cfPot = ColdFusionCheck(pot, coldfusionlimit=coldfusionlimit, coldfusionenergy=coldfusionenergy)
    basinhopping.potential = cfPot
    basinhopping.confCheck.append(cfPot)