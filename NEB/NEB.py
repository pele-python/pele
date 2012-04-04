import numpy as np

def distance_cart(x1, x2):
    return x2 - x1
         
class NEB:
    def __init__(self, initial, final, potential, distance=distance_cart, nimages=10, k=1.0):
        self.distance = distance
        self.potential = potential
        self.k = k
        
        #initialiye coordinate&gradient array
        self.coords = np.zeros([nimages, initial.size])
        self.grad = np.zeros([nimages, initial.size])
        
        #interpolate initial points
        self.interpolate(initial, final, nimages)        
        
        # copy initial and final structure
        self.initial = initial
        self.final = final
        
    
    # do a stupid steepest descent, have to add scipy interface for minimizer first
    def optimize(self, quench=None):
        if(quench==None):
            quench = self.default_quench
        tmp,E = quench(self.coords)
        self.coords = tmp.reshape(self.coords.shape)
        
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

    def getEnergyGradient(self, coords1d):
        coords = coords1d.reshape([10,coords1d.size/10])
        E = 0.
        for i in xrange(0, 10):
            pl = self.initial
            pr = self.final
            if i != 0:
                pl = coords[i - 1, :]
            if i != 9:
                pr = coords[i + 1, :]
            Ei, g = self.NEBForce(coords[i, :], pl, pr)
            self.grad[i,:]=g
            print Ei
            E += Ei
        return E,self.grad.reshape(self.grad.size)
            
    # update force for one image
    def NEBForce(self, p, pl, pr):
            d1 = p - pl
            d2 = pr - p
            t = d1 / np.linalg.norm(d1) + d2 / np.linalg.norm(d2)
            t = t / np.linalg.norm(t)
            
            E, g = self.potential.getEnergyGradient(p)
            
            gpar = g - np.dot(g, t) * t
            gspring = np.dot(self.k * (d2 - d1), t) * t
            print "Energy:", E
            return E, -(gpar + gspring)
    
    # initial interpolation    
    def interpolate(self, initial, final, nimages):
        delta = (final - initial) / (nimages + 1)
        for i in xrange(0, nimages):
            self.coords[i, :] = initial + delta * (i + 1)

            
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
    initial = np.array([1., 2.]) #np.random.random(3)
    final = np.array([2., .7]) #np.random.random(3)
    print "Initial: ", initial
    print "Final: ", final
    #pl.imshow(z)
    
    neb = NEB(initial, final, potential)
    tmp = neb.coords
    
    pl.contourf(x, y, z)
    pl.colorbar()
    pl.plot(tmp[:, 0], tmp[:, 1], 'o-')
    neb.optimize()

        
    tmp = neb.coords
    pl.plot(tmp[:, 0], tmp[:, 1], 'rx-')
    pl.show()
    
        
