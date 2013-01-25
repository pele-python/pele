import numpy as np
from pygmin.transition_states import NEB, NEBPar
from pygmin.transition_states._NEB import distance_cart
from interpolate import InterpolatedPath, interpolate_linear
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
    reinterpolate : integer
        reinterpolate the path to achieve equidistant spacing every so many steps
    adaptive_nimages : bool
        adjust number of images on reinterpolate to match image density
    adaptive_niter : bool
        adjust number of iterations if nimages is adjusted
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
                 k = 100., max_images = -1, image_density=10, iter_density = 10,
                 verbose=0, factor=1., NEBquenchParams=dict(),
                 reinterpolate=0, adaptive_nimages = False, adaptive_niter=False,
                 interpolator=interpolate_linear, distance=distance_cart, parallel=False, ncores=4, **kwargs):
        
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
        self.reinterpolate = reinterpolate
        self._nebclass = NEB
        self._kwargs = kwargs.copy()
        self.k = k
        self.adaptive_images = adaptive_nimages
        self.adaptive_niter = adaptive_niter
        
        if parallel:
            self._kwargs["ncores"]=ncores
            self._nebclass = NEBPar

    @classmethod
    def params(cls, obj = None):
        if obj is None:
            obj = NEBDriver(None, None, None)
        ''' return the the parameters of current instance '''
        params = obj._kwargs.copy()
        params["k"] = obj.k
        params["image_density"] = obj.image_density
        params["max_images"] = obj.max_images
        params["iter_density"] = obj.iter_density
        params["reinterpolate"] = obj.reinterpolate
        params["adaptive_nimages"] = obj.adaptive_images
        params["adaptive_niter"] = obj.adaptive_niter
        
        params["verbose"] = obj.verbose
        params["NEBquenchParams"] = obj.quenchParams.copy()
        # factor is not here, todo, move this out of constuctor
        params["interpolator"] = obj.interpolator
        params["distance"] = obj.distance
        
        if params.has_key("ncores"):
            params["parallel"] = True
         
        return params
    
    def run(self):
                #determine the number of iterations                
        coords1 = self.coords1
        coords2 = self.coords2
            
        path = self.generate_path(coords1, coords2)
        nimages = len(path)

        quenchParams = self.quenchParams.copy()
        
        if quenchParams.has_key("nsteps"):
            niter = quenchParams["nsteps"]
            self.adaptive_niter = False
        else:
            niter = int(self.iter_density * nimages)
            quenchParams["nsteps"] = niter

        #if nimages is already max_images then increasing the number
        #of images with factor will have no effect.  so double the number of steps instead
        if self.factor > 1. and nimages == self.max_images and self.max_images > 0:
            niter *= self.factor
            quenchParams["nsteps"] = niter    
        
        if self.verbose>=0:    
            print "    NEB: nimages", nimages
            print "    NEB: nsteps ", niter
                
        
        if self.reinterpolate > 0:
            quenchParams["nsteps"] = self.reinterpolate    
        
        self.niter = niter
        self.steps_total = 0
        k = self.k
        while True:       
            neb = self._nebclass(path, self.potential, k=k,
                      quenchParams=quenchParams, verbose=self.verbose,
                      distance=self.distance, **self._kwargs)
            self.neb = neb
            #neb.quenchParams["nsteps"]=10
            #print "OPTIMIZING NEB"
            
            neb.events.append(self._process_event)
            res = neb.optimize()
            self.steps_total += res.nsteps
            if res.success or self.steps_total >= self.niter:
                res.nsteps = self.steps_total
                
                if self.verbose >= 0:
                    print "NEB finished after %d steps, rms %e"%(res.nsteps, res.rms)
                return neb
            
            k=neb.k
            distances = []
            for i in xrange(len(path)-1):           
                distances.append(self.distance(res.path[i], res.path[i+1])[0])
            path = self._reinterpolate(res.path, distances)
            
            if self.adaptive_niter:
                self.niter = int(self.iter_density * len(path))
                if self.factor > 1. and len(path) == self.max_images and self.max_images > 0:
                    self.niter *= self.factor    
            if self.verbose >= 1:
                print "NEB reinterpolating path, %d images, niter is %d"%(len(path),self.niter)
        
        
    def generate_path(self, coords1, coords2):
        #determine the number of images to use
        dist, tmp = self.distance(coords1, coords2)
        nimages = int(max(1., dist) * self.image_density * self.factor)
        if self.max_images > 0:
            nimages = min(nimages, self.max_images)
        path = InterpolatedPath(coords1, coords2, nimages, interpolator=self.interpolator)

        return path

    def _reinterpolate(self, path, distances):
        acc_dist = np.sum(distances)
        nimages = len(path)
        if self.adaptive_images:
            nimages = int(max(1., acc_dist) * self.image_density * self.factor)
            
        newpath = []
        newpath.append(path[0].copy())
        
        icur=0
        s_cur = 0.
        s_next=distances[icur]
        
        for i in xrange(1,nimages-1):
            s = float(i)*acc_dist / (nimages-1)
            while s > s_next:
                icur+=1
                s_cur = s_next
                s_next+=distances[icur]
            
            t = (s - s_cur)/(s_next - s_cur)
            newpath.append(self.interpolator(path[icur], path[icur+1], t))
        newpath.append(path[-1].copy())
        return newpath
       
    def _process_event(self, coords=None, energies=None, distances=None, stepnum=None):
        self.update_event(coords=coords, energies=energies,
                       distances=distances, stepnum=stepnum+self.steps_total)