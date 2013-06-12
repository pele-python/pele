from potential import potential

__all__ = ["ColdFusionCheck", "addColdFusionCheck"]

class ColdFusionCheck(potential):
    def __init__(self, potential, coldfusionlimit=-1000, coldfusionenergy=1000000):
        self.potential = potential
        self.coldfusionlimit = coldfusionlimit
        self.coldfusionenergy = coldfusionenergy
        
    def getEnergy(self, coords):
        E = self.potential.getEnergy(coords)
        if(E<self.coldfusionlimit):
            return self.coldfusionenergy
        return E
    
    def getEnergyGradient(self, coords):
        E, grad = self.potential.getEnergyGradient(coords)
        if(E<self.coldfusionlimit):
            grad[:]=0.
        return E, grad
    
    def __call__(self, E, conf, **kwargs):
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