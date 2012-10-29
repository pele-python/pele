import numpy as np
import os.path

import pygmin.defaults as defaults
import pygmin.optimize.quench as quench

__all__ = ["NEB"]

def distance_cart(x1, x2):
    return x2 - x1

class NEB(object):
    """Doubly nudged elastic band implementation

    Parameters
    ----------
    path: iteratable
        iteratable object of all images
    potential :
        the potential object
    distance : callable
        distance function for the elastic band
    k : float
        elastic constant for band
        
    dnep: boolean, optional
        do double nudging, default True
        
    with_springenergy: boolean, optional
        add the spring energy to the total energy of the band, default is False 

    Notes
    -----
    use the Interpolation tools in this package to construct the initial path
    from a starting and ending point
    """
    def __init__(self, path, potential, distance=distance_cart,
                 k=100.0, method="DNEB", with_springenergy=False, dneb=True):
        self.distance = distance
        self.potential = potential
        self.k = k
        nimages = len(path)
        self.nimages = nimages

        self.getEnergyCount = 0
        self.printStateFile = None
        self.iprint = -1


        #initialize coordinate&gradient array
        self.coords = np.zeros([nimages, path[0].size])
        self.energies=np.zeros(nimages)
        self.isclimbing=[]
        for i in xrange(nimages):
            self.isclimbing.append(False)

        # copy initial path
        for i,x in zip(xrange(self.nimages), path):
            self.coords[i,:] = x

        for i in xrange(0,nimages):
            self.energies[i] = potential.getEnergy(self.coords[i,:])
        # the active range of the coords, endpoints are fixed
        self.active = self.coords[1:nimages-1,:]
        
        self.dneb = dneb
        self.with_springenergy = with_springenergy

    def optimize(self, quenchRoutine=None,
                 **kwargs):
        """
        Optimize the band

        Note: the potential for the NEB optimization is not Hamiltonian.  This
            means that there is no meaningful energy associated with the
            potential.  Therefore, during the optimization, we can do gradient
            following, but we can't rely on the energy for, e.g. determining
            step size, which is the default behavior for many optimizers.  This
            can be worked around by choosing a small step size and a large
            maxErise, or by using an optimizer that uses only gradients.

            scipy.lbfgs_b seems to work with NEB pretty well, but lbfgs_py and
            mylbfgs tend to fail.  If you must use one of those try, e.g.
            maxErise = 1., maxstep=0.01, tol=1e-2

        :quenchRoutine: quench algorithm to use for optimization.

        :quenchParams: parameters for the quench """ 
        if quenchRoutine is None:
            quenchRoutine = defaults.NEBquenchRoutine 
        #combine default and passed params.  passed params will overwrite default 
        quenchParams = dict(defaults.NEBquenchParams.items() +
                            kwargs.items() )

        if quenchParams.has_key("iprint"):
            self.iprint = quenchParams["iprint"]

        tmp,E,rms,tmp4 = quenchRoutine(
                    self.active.reshape(self.active.size), self.getEnergyGradient,
                    **quenchParams)
        print "neb rms", rms
        self.active[:,:] = tmp.reshape(self.active.shape)
        for i in xrange(0,self.nimages):
            self.energies[i] = self.potential.getEnergy(self.coords[i,:])

    def getEnergyGradient(self, coords1d):
        """
        Calculates the gradient for the whole NEB. only use force based minimizer!

        coords1d:
            coordinates of the whole neb active images (no end points)
        """
        # make array access a bit simpler, create array which contains end images
        tmp = self.coords.copy()
        tmp[1:self.nimages-1,:] = coords1d.reshape(self.active.shape)
        grad = np.zeros(self.active.shape)

        # calculate real energy and gradient along the band. energy is needed for tangent
        # construction
        realgrad = np.zeros(tmp.shape)
        for i in xrange(1, self.nimages-1):
            self.energies[i], realgrad[i,:] = self.potential.getEnergyGradient(tmp[i,:])

        # the total energy of images, band is neglected
        E = sum(self.energies)
        Eneb = 0
        # build forces for all images
        for i in xrange(1, self.nimages-1):
            En, grad[i-1,:] = self.NEBForce(
                    self.isclimbing[i],
                    [self.energies[i],tmp[i, :]],
                    [self.energies[i-1],tmp[i-1, :]],
                    [self.energies[i+1],tmp[i+1, :]],
                    realgrad[i,:]
                    )
            Eneb += En
        if self.iprint > 0:
            if self.getEnergyCount % self.iprint == 0:
                self.printState()
        self.getEnergyCount += 1
        #print "ENeb = ", Eneb
        return E+Eneb, grad.reshape(grad.size)
        #return 0., grad.reshape(grad.size)

    def tangent_old(self, central, left, right):
        """
        Old tangent construction based on average of neighbouring images

        coords1d:
            coordinates of the whole neb active images (no end points)
        """
        d1 = self.distance(central[1], left[1])
        d2 = self.distance(right[1], central[1])
        t = d1 / np.linalg.norm(d1) + d2 / np.linalg.norm(d2)
        return t / np.linalg.norm(t)

    def tangent(self, central, left, right):
        """
        New uphill tangent formulation

        The method was  described in
        "Improved tangent estimate in the nudged elastic band method for finding
        minimum energy paths and saddle points"
        Graeme Henkelman and Hannes Jonsson
        J. Chem. Phys 113 (22), 9978 (2000)

        central:
            central image energy and coordinates [E, coords]
        left:
            left image energy and coordinates [E, coords]
        right:
            right image energy and coordinates [E, coords]
        """
        tleft = self.distance(central[1], left[1])
        tright = self.distance(right[1],  central[1])
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

    def NEBForce(self, isclimbing, image, left, right, greal):
        """
        Calculate NEB force for 1 image. That contains projected real force and spring force.

        The current implementation is the DNEB (doubly nudged elastic band) as described in

        "A doubly nudged elastic band method for finding transition states"
        Semen A. Trygubenko and David J. Wales
        J. Chem. Phys. 120, 2082 (2004); doi: 10.1063/1.1636455

        """
        if(isclimbing):
            return greal - 2.*np.dot(greal, t) * t

        # construct tangent vector
        p = image[1]
        pl = left[1]
        pr = right[1]

        #d1 = self.distance(image[1], left[1])
        #d2 = self.distance(right[1], image[1])


        t = self.tangent(image,left,right)
        if True:
            import _NEB_utils
            E, g_tot = _NEB_utils.neb_force(t,greal, self.k, self.dneb)
            if self.with_springenergy:
                return E, g_tot
            else:
                return 0., g_tot
        else:
    
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
                                    
            g_tot = gperp + gs_par
    
            if(self.dneb):
                # double nudging
                g_tot += gs_perp - np.dot(gs_perp,gperp)*gperp/np.dot(gperp,gperp)
            
            if(self.with_springenergy):
                E = 0.5 / self.k * np.dot(gspring, gspring)
            else:
                E = 0.
    
            return E, g_tot

    def MakeHighestImageClimbing(self):
        """
        Make the image with the highest energy a climbing image
        """
        emax = max(self.energies)
        for i in xrange(1,len(self.energies)-1):
            if(abs(self.energies[i]-emax)<1e-10):
                self.isclimbing[i] = True

    def MakeAllMaximaClimbing(self):
        """
        Make all maxima along the neb climbing images.
        """
        for i in xrange(1,len(self.energies)-1):
            if(self.energies[i] > self.energies[i-1] and self.energies[i] > self.energies[i+1]):
                self.isclimbing[i] = True

    def printState(self):
        """
        print the current state of the NEB.  Useful for bug testing
        """
        if self.printStateFile is None:
            #get a unique file name
            for n in range(500):
                self.printStateFile = "neb.EofS.%04d" % (n)
                if not os.path.isfile(self.printStateFile):
                    break
            print "NEB energies will be written to", self.printStateFile
        #fname = "neb.EofS.%04d" % (self.getEnergyCount)
        #fname = "neb.EofS.all"
        with open(self.printStateFile, "a") as fout:
            fout.write( "\n\n" )
            fout.write("#neb step: %d\n" % (self.getEnergyCount))
            S = 0.
            for i in range(self.nimages):
                if i == 0:
                    dist = 0.
                else:
                    dist = self.distance( self.coords[i,:], self.coords[i-1,:] )
                    dist = np.linalg.norm(dist)
                S += dist
                #print "S",S, "E",self.energies[i], "dist", dist
                fout.write("%f %g\n" % (S, self.energies[i]))





import nebtesting as test

if __name__ == "__main__":
    import pylab as pl
    from interpolate import InterpolatedPath
    from pygmin import defaults
    from pygmin.optimize import quench
    defaults.NEBquenchRoutine = quench.lbfgs_py
    defaults.NEBquenchParams["iprint"]=1
    defaults.NEBquenchParams["debug"]=True
    defaults.NEBquenchParams["maxErise"]=0.1
    
    x = np.arange(.5, 5., .05)
    y = np.arange(.5, 5., .05)
    z = np.zeros([len(x), len(y)])
    potential = test.leps()
    for i in range(0, len(x)):
        for j in range(0, len(y)):
                z[j, i] = potential.getEnergy([x[i], y[j]])
    print "done"
    #z[z>0.] = 0.
    #pl.imshow(z)
    #pl.show()
    initial = np.array([.75, 2.]) #np.random.random(3)
    final = np.array([2., .75]) #np.random.random(3)
#    from pygmin.optimize import quench
#    print "quench initial"
#    ret = quench.lbfgs_py(initial, potential.getEnergyGradient)
#    initial = ret[0]
#    print "quench final"
#    ret = quench.quench(final, potential.getEnergyGradient)
#    final = ret[0]
#    print "done with quenching"
#    print initial, final
    #print "Initial: ", initial
    #print "Final: ", final
    #pl.imshow(z)

    neb = NEB(InterpolatedPath(initial, final, 20) ,potential, k=1000)
    tmp = neb.coords
    energies_interpolate = neb.energies.copy()
    pl.figure()
    pl.subplot(1,2,1)
    pl.subplots_adjust(wspace=0.3, left=0.05, right=0.95, bottom = 0.14)

    pl.title("path")
    #pl.contourf(x, y, z)
    pl.pcolor(x, y, z, vmax=-0.5, cmap=pl.cm.PuBu)
    pl.colorbar()
    pl.plot(tmp[:, 0], tmp[:, 1], 'ko-')
    print "optimizing NEB"
    neb.optimize()#quenchRoutine=quench.fire)
    print "done"
    tmp = neb.coords
    pl.plot(tmp[:, 0], tmp[:, 1], 'ro-')
    pl.xlabel("x")
    pl.ylabel("y")
    pl.axis(xmin=0.5, xmax=2.5, ymin=0.5, ymax=2.5)

    pl.subplot(1,2,2)
    pl.title("energy")
    pl.plot(energies_interpolate, 'ko-', label="interpolate")
    pl.plot(neb.energies, 'ro-', label="neb")
    pl.xlabel("image")
    pl.ylabel("energy")
    pl.legend(loc='best')
    pl.show()
