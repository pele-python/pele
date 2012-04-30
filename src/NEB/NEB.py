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
    def __init__(self, initial, final, potential, distance=distance_cart, nimages=20, k=100.0):
        self.distance = distance
        self.potential = potential
        self.k = k
        self.nimages = nimages
        
        #initialiye coordinate&gradient array
        self.coords = np.zeros([nimages, initial.size])
        self.grad = np.zeros([nimages-2, initial.size])
        self.energies=np.zeros(nimages)
        self.isclimbing=[]
        for i in xrange(nimages):
            self.isclimbing.append(False)
            
        #interpolate initial points
        self.interpolate(initial, final, nimages)        
        
        # copy initial and final structure
        self.coords[0,:] = initial
        self.coords[-1,:] = final
        for i in xrange(0,nimages):
            self.energies[i] = potential.getEnergy(self.coords[i,:])
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
        for i in xrange(0,self.nimages):
            self.energies[i] = self.potential.getEnergy(self.coords[i,:])
        
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

    def default_quench_fire(self, coords):
        import optimize.fire as fire
        tmp = coords.reshape(coords.size)
        opt = fire.Fire(tmp, self.getEnergyGradient,dtmax=0.1, dt=0.01, maxmove=0.01)
        opt.run()
        return opt.coords, 0.0
    
    # Calculate gradient for the while NEB
    def getEnergyGradient(self, coords1d):
        # make array access a bit simpler
        tmp = self.coords.copy()
        tmp[1:self.nimages-1,:] = coords1d.reshape(self.active.shape)
        grad = self.grad.copy()
        
        # calculate real energy and gradient along the band
        self.realgrad = np.zeros(tmp.shape)
        for i in xrange(1, self.nimages-1):
            self.energies[i], self.realgrad[i,:] = self.potential.getEnergyGradient(tmp[i,:])

        E = sum(self.energies)
        
        # build forces for all images
        for i in xrange(1, self.nimages-1):
            g = self.NEBForce(
                    [self.energies[i],tmp[i, :]],
                    [self.energies[i-1],tmp[i-1, :]],
                    [self.energies[i+1],tmp[i+1, :]],
                    self.realgrad[i,:],
                    self.isclimbing[i])
            
            grad[i-1,:] = g
        return E,grad.reshape(self.grad.size)
            
    # old average tangent formulation
    def tangentaa(self, central, left, right):
        d1 = central[1] - left[1]
        d2 = right[1] - central[1]
        t = d1 / np.linalg.norm(d1) + d2 / np.linalg.norm(d2)
        return t / np.linalg.norm(t)
        
    # new uphill tangent formulation
    def tangent(self, central, left, right):
        tleft = (central[1] - left[1])        
        tright = (right[1] - central[1])
        vmax = max(abs(central[0] - left[0]), abs(central[0] - right[0]))
        vmin = max(abs(central[0] - left[0]), abs(central[0] - right[0]))
        
        # special interpolation treatment for maxima/minima
        if (central[0] >= left[0] and central[0] >= right[0]) or (central[0] <= left[0] and central[0] <= right[0]):
            if(left[0] > right[0]):
                t = vmax * tleft + vmin*tright
            else:
                t = vmin * tleft + vmax*tright
        # left is higher, take this one
        elif (left[0] > right[0]):
            t = tleft
        # otherwise take right
        else: 
            t = tright

        return t / np.linalg.norm(t)            
    
    # update force for one image
    def NEBForce(self, image, left, right, greal, isclimbing):
            # construct tangent vector, TODO: implement newer method
            p = image[1]
            pl = left[1]
            pr = right[1]
            
            d1 = image[1] - left[1]
            d2 = right[1] - image[1]
        
            
            t = self.tangent(image,left,right)
            
            # project out parallel part
            gperp = greal - np.dot(greal, t) * t
            # calculate parallel spring force and energy
            #gspring = -self.k * (np.linalg.norm(d2) - np.linalg.norm(d1)) * t
            # this is the spring
            gspring = -self.k*(pl + pr - 2.*p)
            # the parallel part
            gs_par = np.dot(gspring,t)*t
            # perpendicular part
            gs_perp = gspring - gs_par
            # double nudging
            gstar = gs_perp - np.dot(gs_perp,gperp)*gperp/np.dot(gperp,gperp)            
            
            if(isclimbing):
                return greal - 2.*np.dot(greal, t) * t
            return (gperp + gs_par + gstar)
    
    # initial interpolation    
    def interpolate(self, initial, final, nimages):
        delta = (final - initial) / (nimages-1)
        for i in xrange(1, nimages):
            self.coords[i, :] =  initial + delta * i
            
    def MakeClimbingImage(self):
        emax = max(self.energies)
        for i in xrange(0,len(self.energies)):
            if(abs(self.energies[i]-emax)<1e-10):
                self.isclimbing[i] = True
            
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
    #print "Initial: ", initial
    #print "Final: ", final
    #pl.imshow(z)
    
    neb = NEB(initial, final, potential, nimages=20, k=1000)
    tmp = neb.coords
    
    pl.contourf(x, y, z)
    pl.colorbar()
    pl.plot(tmp[:, 0], tmp[:, 1], 'o-')
    neb.optimize()

    tmp = neb.coords
    pl.plot(tmp[:, 0], tmp[:, 1], 'ro-')
    pl.show()
    print "bla"
    pl.plot(neb.energies)
    pl.show()
    print "bla2" 
    
        
