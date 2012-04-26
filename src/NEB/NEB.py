import numpy as np

def distance_cart(x1, x2):
    return x2 - x1

class NEB:
    """Nudged elastic band implementation

    intial: 
        first point of band
        
    final: last point of band
    
    distance  (distance_cart):
        distance function for the elastic band
    
    nimages (10):
        number of moving images for the band. The number includes 
        the endpoints and must be bigger than 2
        
    k (1.0):
        elastic constant for band
        
    """
    def __init__(self, initial, final, potential, distance=distance_cart, nimages=10, k=20.0):
        self.distance = distance
        self.potential = potential
        self.k = k
        self.nimages = nimages
        
        #initialiye coordinate&gradient array
        self.coords = np.zeros([nimages, initial.size])
        self.grad = np.zeros([nimages-2, initial.size])
        
        #interpolate initial points
        self.interpolate(initial, final, nimages)        
        
        # copy initial and final structure
        self.coords[0,:] = initial
        self.coords[-1,:] = final
        
        # the active range of the coords, endpoints are fixed
        self.active = self.coords[1:nimages-1,:]  
    
    # do a stupid steepest descent, have to add scipy interface for minimizer first
    """
        Optimize the band
    
        quench (None):
            quench algorithm to use for optimization. If None is given,
            default_quench is used,
    """
    def optimize(self, quench=None):
        if(quench==None):
            quench = self.default_quench
        tmp,E = quench(self.active)
        self.active[:,:] = tmp.reshape(self.active.shape)
        
    def default_quench(self, coords):
        import scipy.optimize.lbfgsb
        #newcoords, newE = steepest_descent.steepestDescent(potential.getEnergyGradient, coords, 100)
        tmp = coords.reshape(coords.size)
        newcoords, newE, dictionary = scipy.optimize.fmin_l_bfgs_b(self.getEnergyGradient, tmp, iprint=-1, pgtol=1e-3)

        warnflag = dictionary['warnflag']
        if warnflag > 0:
            print "warning: problem with quench: ",
            if warnflag == 1:
                print "too many function evaluations"
            else:
                print dictionary['task']            
        return newcoords, newE

    # Calculate gradient for the while NEB
    def getEnergyGradient(self, coords1d):
        # make array access a bit simpler
        tmp = self.coords.copy()
        tmp[1:self.nimages-1,:] = coords1d.reshape(self.active.shape)
        grad = self.grad.copy()
        E = 0.
        # build forces for all images
        for i in xrange(1, self.nimages-1):
            Ei, g = self.NEBForce(tmp[i, :], tmp[i-1, :], tmp[i+1, :])
            grad[i-1,:] = g
            E += Ei
        return E,grad.reshape(self.grad.size)
            
    # update force for one image
    def NEBForce(self, p, pl, pr):
            # construct tangent vector, TODO: implement newer method
            d1 = p - pl
            d2 = pr - p
            t = d1 / np.linalg.norm(d1) + d2 / np.linalg.norm(d2)
            t = t / np.linalg.norm(t)
            
            # get real gradient for image
            E, g = self.potential.getEnergyGradient(p)
            
            # project out parallel part
            gpar = g - np.dot(g, t) * t
            # calculate parallel spring force and energy
            gspring = self.k * (np.linalg.norm(d2) - np.linalg.norm(d1)) * t
            E+=0.5*self.k*(np.linalg.norm(d1)-np.linalg.norm(d2))**2
            
            return E, -(gpar + gspring)
    
    # initial interpolation    
    def interpolate(self, initial, final, nimages):
        delta = (final - initial) / (nimages-1)
        for i in xrange(1, nimages):
            self.coords[i, :] =  initial + delta * i

            
import nebtesting as test

if __name__ == "__main__":
    import pylab as pl
    x = np.arange(.5, 3., .1)
    y = np.arange(.5, 3., .1)
    z = np.zeros([len(x), len(y)])
    potential = test.leps()
    for i in range(0, len(x)):
        for j in range(0, len(y)):
                z[j, i] = potential.getEnergy([x[i], y[j]])
    print "done"
    #pl.imshow(z)
    #pl.show()
    initial = np.array([.75, 2.]) #np.random.random(3)
    final = np.array([2., .75]) #np.random.random(3)
    print "Initial: ", initial
    print "Final: ", final
    #pl.imshow(z)
    
    neb = NEB(initial, final, potential, nimages=20)
    tmp = neb.coords
    
    pl.contourf(x, y, z)
    pl.colorbar()
    pl.plot(tmp[:, 0], tmp[:, 1], 'o-')
    neb.optimize()

        
    tmp = neb.coords
    pl.plot(tmp[:, 0], tmp[:, 1], 'ro-')
    pl.show()
    
        
