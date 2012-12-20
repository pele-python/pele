from pygmin.optimize import lbfgs_scipy as quench_lbfgs


class Minimum:
    #use this to define a minimum, 
    #and to how to determine whether two minima are the same
    def __init__(self, energy, coords):
        self.energy = energy
        #self.coords = coords.copy()
    def isEqual(self, min2):
        return abs(self.energy - min2.energy) < 1e-4 


class DoQuenching:
    def __init__(self, potential, coords, quench = quench_lbfgs):
        self.potential = potential
        self.quench = quench
        
        #do initial quench
        ret = self.quench(coords, potential.getEnergyGradient)
        energy = ret[1]
        self.minimum = Minimum(energy, coords)
        
        #save quenched energy
        
        self.nrejected = 0
        self.ntot = 0
    
    def acceptReject(self, oldE, newE, old_coords, new_coords):
        #do quenching
        coords = new_coords
        ret = self.quench(coords, self.potential.getEnergyGradient)
        energy = ret[1]
        newmin = Minimum(energy, coords)
        
        accepted = self.minimum.isEqual( newmin )
        self.ntot += 1
        if not accepted: 
            self.nrejected += 1
            print "step quenched to new minimum", self.minimum.energy, newmin.energy, accepted
        

        # here do whatever else you want.
        # e.g. quench with different quenchers, etc.
                
        return accepted
