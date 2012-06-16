from potential import potential

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
            print "Cold fusion detected"
            grad[:]=0.
            E = self.coldfusionenergy
        return E, grad