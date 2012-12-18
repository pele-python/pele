import numpy as np
from pygmin.potentials import BasePotential


class XYModel(BasePotential):
    """
    1d xymodel
    """
    def __init__(self, nspins, phi=1.0, phases=None):
        self.nspins = nspins
        
        if phases == None:
            self.phases = np.random.uniform(-phi,phi, self.nspins)
        else:
            self.phases = phases
        
        self.periodic = True
        
    def getEnergy(self, angles):
        #do internal energies first
        E = np.sum( np.cos( -angles[0:-1] + angles[1:] + self.phases[0:-1] ) )
        
        #now do boundary condition
        if self.periodic:
            E += np.cos( -angles[-1] + angles[0] + self.phases[-1] )
        
        E = self.nspins - E #/ self.nspins
        return E

    def getEnergyGradient(self, angles):
        from numpy import cos, sin
        
        grad = np.zeros(self.nspins)
        
        anglediff = np.zeros(self.nspins)
        anglediff[0:-1] = -angles[0:-1] + angles[1:] + self.phases[0:-1]
        if self.periodic:
            anglediff[-1] = -angles[-1] + angles[0] + self.phases[-1]
        
        #first do interior spins
        grad[1:-1]= sin( anglediff[0:-2] ) - sin( anglediff[1:-1] ) 
        
        #now do spins on the boundaries
        grad[0] =  sin( anglediff[-1] ) - sin( anglediff[0] )
        grad[-1] = sin( anglediff[-2] ) - sin( anglediff[-1] )
        
        E = np.sum( cos( anglediff[0:-1] ) )
        E += cos( anglediff[-1] )
        E = self.nspins - E# / self.nspins



        #grad /= self.nspins
        return E, grad

#    def getEnergyGradient(self, angles):
#        e = self.getEnergy(angles)
#        grad = self.getGradient(angles)
#        return e, grad


def test_basin_hopping(pot, angles):
    from pygmin.basinhopping import BasinHopping
    from pygmin.takestep.displace import RandomDisplacement
    from pygmin.takestep.adaptive import AdaptiveStepsize
    
    takestep = RandomDisplacement(stepsize = np.pi/4)
    takestepa = AdaptiveStepsize(takestep, frequency = 20)
    
    bh = BasinHopping( angles, pot, takestepa, temperature = 1.01)
    bh.run(400)

def test():
    pi = np.pi
    nspins = 6
    
    #phases = np.zeros(nspins)
    pot = XYModel(nspins, phi = np.pi/8.) #, phases=phases)
        
    angles = np.random.uniform(-pi/4, pi/4, nspins)
    print angles

    e = pot.getEnergy(angles)
    print e
    
    print "numerical gradient"
    ret = pot.getEnergyGradientNumerical(angles)
    print ret[1]
    print "analytical gradient"
    ret2 = pot.getEnergyGradient(angles)
    print ret2[1]
    print ret[0]
    print ret2[0]
    
    
    #try a quench
    from pygmin.optimize.quench import mylbfgs
    ret = mylbfgs(angles, pot.getEnergyGradient)
    
    print "quenched e = ", ret[1]
    print ret[0]
    
    test_basin_hopping(pot, angles)

if __name__ == "__main__":
    test()