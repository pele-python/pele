import numpy as np
from pygmin.transition_states import NEB, NEBPar
from pygmin.transition_states._NEB import distance_cart
from interpolate import InterpolatedPath
from pygmin.utils.events import Signal

#def calc_neb_dist(coords, nimages, dist=True, grad=False):
#    d_left = np.zeros(coords.shape)
#    coords = coords.reshape([-1,nimages])
#    for i in xrange(nimages):
        

#def create_NEB(pot, coords1, coords2, image_density=10, max_images=40,
#                iter_density=15, 
#                NEBquenchParams=dict(),
#                interpolator=None,
#                verbose=False, factor=1, parallel=False, ncores=4, **NEBparams):

class NEBDriver(object):
    ''' driver class for NEB
    
    The NEBDriver wraps calls for NEB from LocalConnect. The driver class is responsible for setting
    up the initial interpolation and optimizing the band.
    
    Parameters:
    -----------
    potential :
        the potential object

    coords1, coords2 : array
        the structures to connect with the band
        
    max_images : int
        the maximum number of NEB images
    image_density : float
        how many NEB images per unit distance to use.
    iter_density : float
        how many optimization iterations per unit distance to use.
    factor : int
        The number of images is multiplied by this factor.  If the number of 
        images is already at it's maximum, then the number of iterations is 
        multiplied by this factor instead
    verbose : integer
    parallel : bool
        if True, then use class NEBPar to evaluate the image potentials in parallel
    ncores : int
        the number of cores to use.  Ignored if parallel is False
    interpolator : callable, optional
        the function used to do the path interpolation for the NEB
        
    See Also
    ---------
    NEB
    InterpolatedPath
        
    '''
    
    def __init__(self, potential, coords1, coords2,
                 max_images = -1, image_density=10, iter_density = 10,
                 verbose=-1, factor=1., NEBquenchParams=dict(),
                 interpolator=None, distance=distance_cart, parallel=False, ncores=4, **kwargs):
        
        self.potential = potential
        self.interpolator = interpolator
        
        self.verbose = verbose
        self.distance = distance
        self.factor = factor
        self.max_images = max_images
        self.image_density = image_density
        self.iter_density = iter_density
        self.update_event = Signal()
        self.quenchParams=NEBquenchParams
        self.coords1 = coords1
        self.coords2 = coords2   
             
        self._nebclass = NEB
        self._kwargs = kwargs.copy()
        
        if parallel:
            self._kwargs["ncores"]=ncores
            self._nebclass = NEBPar

    def run(self):
                #determine the number of iterations                
        coords1 = self.coords1
        coords2 = self.coords2
            
        path = self.generate_path(coords1, coords2)
        nimages = len(path)

        quenchParams = self.quenchParams.copy()
        
        if quenchParams.has_key("nsteps"):
            niter = quenchParams["nsteps"]
        else:
            niter = int(self.iter_density * nimages)
            quenchParams["nsteps"] = niter

        #if nimages is already max_images then increasing the number
        #of images with factor will have no effect.  so double the number of steps instead
        if self.factor > 1. and nimages == self.max_images and self.max_images > 0:
            niter *= self.factor
            quenchParams["nsteps"] = niter    
        
        if self.verbose:    
            print "    NEB: nimages", nimages
            print "    NEB: nsteps ", niter
                
        neb = self._nebclass(path, self.potential,
                  quenchParams=quenchParams, verbose=self.verbose,
                  distance=self.distance, **self._kwargs)
        
        #neb.quenchParams["nsteps"]=10
        #print "OPTIMIZING NEB"
        neb.events.append(self.update_event)
        neb.optimize()
        
        return neb        
       
    def generate_path(self, coords1, coords2):
        #determine the number of images to use
        dist, tmp = self.distance(coords1, coords2)
        nimages = int(max(1., dist) * self.image_density * self.factor)
        if self.max_images > 0:
            nimages = min(nimages, self.max_images)
        path = InterpolatedPath(coords1, coords2, nimages, interpolator=self.interpolator)

        return path
          