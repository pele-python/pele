def interpolate_linear(initial, final, t):
    '''
        Linear interpolation
        
        Parameters
        ----------
        initial: 
            first point
        final:
            second point
        t: float
            weight for interpolation [0.,1.]
    '''
    return (1.-t)*initial + t*final

def InterpolatedPathDensity(initial, final, distance, density=10., **kwargs):
    """
    Return a InterpolatedPath object with the appropriate 
    number of images corresponding to the given density
    
    Params
    ------
    density : 
        number of images per "unit" separation
    kwargs : 
        parameters passed to InterpolatedPath
    """
    if distance > 1:
        nimages = int(density * distance)
    else:
        nimages = int(density)
    return InterpolatedPath(initial, final, nimages, **kwargs)
    

class InterpolatedPath(object):
    '''
        Wraps interpolation that it can be accessed like a list / array without storing the nodes
        
        Parameters
        ----------
        intitial:
            first point
            
        final:
            second point
            
        nimages: integer
            number of images
            
        interpolator: callable
            interpolation algorithm
            
        Examples
        --------
        
        >>> import numpy as np
        >>>
        >>> path = InterpolatedPath(np.array([0., 1.]), np.array([1.,0.]), 10)
        >>> print len(path)  # prints 10
        >>> print path[5]    # access 5th frame
        >>> 
        >>> for x in path: # it is also iteratable
        >>>     print x
    '''
    initial=None
    final=None
    nimages=0
    interpolator=None
    current = 0
    
    def __init__(self, initial, final, nimages, interpolator=interpolate_linear):
        self.initial=initial
        self.final=final
        self.nimages = nimages
        self.interpolator = interpolator
        
    
    def __len__(self):
        return self.nimages
    
    def __getitem__(self, index):
        assert index>=0 and index<=self.nimages
        return self.interpolator(self.initial, self.final, float(index) / float(self.nimages-1))           
         
    #required iterable elements
    class Iterator(object):
        def __init__(self, path):
            self.path = path
            self.index = -1
            
        def __iter__(self):
            return self
        
        def next(self):
            if self.index == self.path.nimages-1:
                raise StopIteration
            self.index += 1
            return self.path.__getitem__(self.index)
    
    def __iter__(self):
        return self.Iterator(self)
    
if __name__ == "__main__":
    path = InterpolatedPath(0., 1., 10)
    print len(path)
    print path[5]
    
    i=0
    for x in path:
        print i,x
        i+=1