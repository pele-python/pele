import tempfile

from pygmin.landscape import DoubleEndedConnect, DoubleEndedConnectPar
from pygmin import basinhopping
from pygmin.storage import Database
from pygmin.takestep import RandomDisplacement, AdaptiveStepsizeTemperature
from pygmin.utils.xyz import write_xyz
from pygmin.optimize import mylbfgs

__all__ = ["BaseParameters", "Parameters", "dict_copy_update", "BaseSystem"]

#class NotImplemented(BaseException):
#    """
#    The exception to return if there is a feature
#    in the System is not implemented yet
#    """
#    pass


class BaseParameters(dict):
    """define a dictionary who's values can be accessed like attributes
    
    if `params` is a BaseParameters object then this command::
    
        base_params["key"] = value
    
    is the same as::
    
        base_params.key = value
        
    This only works for keys that are strings
    
    """
    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name)

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

class Parameters(BaseParameters):
    """Define the parameter tree for use with BaseSystem class"""
    def __init__(self):
        self["database"] = BaseParameters()
        self["basinhopping"] = BaseParameters()
        self["takestep"] = BaseParameters()
        self.structural_quench_params = BaseParameters()
        self.gui = BaseParameters()
        
        
        self.double_ended_connect = BaseParameters()
        self.double_ended_connect.local_connect_params = BaseParameters()

        self.double_ended_connect.local_connect_params.pushoff_params = BaseParameters()
        self.double_ended_connect.local_connect_params.pushoff_params.quenchParams = BaseParameters()

        self.double_ended_connect.local_connect_params.tsSearchParams = BaseParameters()
        self.double_ended_connect.local_connect_params.NEBparams = BaseParameters()
        self.double_ended_connect.local_connect_params.NEBparams.NEBquenchParams = BaseParameters()
        
        self.double_ended_connect.local_connect_params.tsSearchParams.lowestEigenvectorQuenchParams = BaseParameters()
        self.double_ended_connect.local_connect_params.tsSearchParams.tangentSpaceQuenchParams = BaseParameters()     


def dict_copy_update(dict1, dict2):
    """return a new dictionary from the union of dict1 and dict2.  if there
    are conflicts, take the value in dict2"""
    newdict = dict1.copy()
    newdict.update(dict2)
    return newdict

class BaseSystem(object):
    """
    this class defines a base class for a System object
        
    Notes
    -------------
    The following functions need to be overloaded for running
    the given routines

    
    Global Optimization::

    1. get_potential : required
    #. get_takestep : optional
    #. get_random_configuration : optional
    #. get_compare_exact : optional

    Connecting Minima and Transition State Searches::

    1. get_potential : required
    #. get_mindist : required
    #. get_orthogonalize_to_zero_eigenvectors : required
    #. get_compare_exact : optional, recommended
    
    GUI::
        
    1. all of the above functions are required, plus
    #. draw : required
    #. smooth_path : required
    #. load_coords_pymol : recommended

    additionally, it's a very good idea to specify the accuracy in the 
    database using self.params.database.accuracy
    
    See the method documentation for more information and relevant links
    
    """
    def __init__(self, *args, **kwargs):
        self.params = Parameters()
        
        tsSearchParams = self.params.double_ended_connect.local_connect_params.tsSearchParams
        tsSearchParams.lowestEigenvectorQuenchParams.nsteps = 100
        
#        self.params.double_ended_connect.local_connect_params.NEBparams.NEBquenchParams.maxErise = 1e50

    def get_potential(self):
        """return the potential object
        
        See Also
        --------
        pygmin.potentials
        """
        raise NotImplementedError

    def get_random_configuration(self):
        """a starting point for basinhopping, etc."""
        raise NotImplementedError
    
    def get_random_minimized_configuration(self):
        coords = self.get_random_configuration()
        quencher = self.get_minimizer()
        return quencher(coords)
    
    def get_minimizer(self):
        """return a function to minimize the structure"""
        pot = self.get_potential()
        params = self.params.structural_quench_params
        return lambda coords: mylbfgs(coords, pot.getEnergyGradient, **params)
    
    def get_compare_exact(self):
        """object that returns True if two structures are exact.
        
            true_false = compare_exact(min1, min2)
        
        See Also
        --------
        pygmin.mindist
        """
        raise NotImplementedError

    def get_compare_minima(self):
        """a wrapper for compare exact so in input can be in 
        Minimum Form"""
        compare = self.get_compare_exact()
        if compare is None:
            return None
        return lambda m1, m2: compare(m1.coords, m2.coords)
    
    def create_database(self, *args, **kwargs):
        """return a new database object
        
        See Also
        --------
        pygmin.storage
        """
        kwargs = dict_copy_update(self.params["database"], kwargs)        
        #note this syntax is quite ugly, but we would like to be able to 
        #create a new database by passing the filename as the first arg, 
        #not as a kwarg.  
        if len(args) > 1:
            raise ValueError("create_database can only take one non-keyword argument")
        if len(args) == 1:
            if "db" not in kwargs:
                kwargs["db"] = args[0]

        #get a routine to compare the minima as exact
        try:
            if not "compareMinima" in kwargs:
                try:
                    compare_minima = self.get_compare_minima()
                    kwargs["compareMinima"] = compare_minima
                except NotImplementedError:
                    pass
        except NotImplementedError:
            #compareMinima is optional
            pass

        return Database(**kwargs)
    

    def get_takestep(self, stepsize=0.6, **kwargs):
        """return the takestep object for use in basinhopping, etc.
        
        default is random displacement with adaptive step size 
        adaptive temperature
        
        See Also
        --------
        pygmin.takestep
        """
        kwargs = dict_copy_update(self.params["takestep"], kwargs)        
        takeStep = RandomDisplacement(stepsize=stepsize)
        tsAdaptive = AdaptiveStepsizeTemperature(takeStep, **kwargs)
        return tsAdaptive

    def get_basinhopping(self, database=None, takestep=None, coords=None, add_minimum=None, **kwargs):
        """return the basinhopping object with takestep
        and accept step already implemented
        
        See Also
        --------
        pygmin.basinhopping
        """
        kwargs = dict_copy_update(self.params["basinhopping"], kwargs)
        pot = self.get_potential()
        if coords is None:
            coords = self.get_random_configuration()
        if takestep is None:
            takestep = self.get_takestep()
        if add_minimum is None:
            if database is None:
                database = self.create_database()
            add_minimum = database.minimum_adder()
        bh = basinhopping.BasinHopping(coords, pot, takestep, storage=add_minimum, **kwargs)
        return bh

    def get_mindist(self):
        """return a mindist object that is callable with the form
        
        dist, X1new, X2new = mindist(X1, X2)
        
        Notes
        -----
        the mindist object returns returns the best alignment between two
        configurations, taking into account all global symmetries
        
        See Also
        --------
        pygmin.mindist
         
        """
        raise NotImplementedError
    
    def get_orthogonalize_to_zero_eigenvectors(self):
        """return a which makes a vector orthogonal to the known zero
        eigenvectors (the eigenvectors with zero eigenvalues.  It should
        be callable with the form
        
            vec = orthogVec(vec, coords)
        
        See Also
        --------
        pygmin.transition_states
        """
        raise NotImplementedError
    
    def get_double_ended_connect(self, min1, min2, database, parallel=False, **kwargs):
        """return a DoubleEndedConnect object
    
        See Also
        --------
        pygmin.landscape
        """
        kwargs = dict_copy_update(self.params["double_ended_connect"], kwargs)
        pot = self.get_potential()
        mindist = self.get_mindist()
        
        #attach the function which orthogonalizes to known zero eigenvectors.
        #This is amazingly ugly
        try:
            kwargs["local_connect_params"]["tsSearchParams"]["orthogZeroEigs"]
        except KeyError:
            if "local_connect_params" in kwargs:
                lcp = kwargs["local_connect_params"]
            else:
                lcp = kwargs["local_connect_params"] = BaseParameters()
            
            if "tsSearchParams" in lcp:
                tssp = lcp["tsSearchParams"]
            else:
                tssp = lcp["tsSearchParams"] = BaseParameters()
            
            if not "orthogZeroEigs" in tssp:
                tssp["orthogZeroEigs"] = self.get_orthogonalize_to_zero_eigenvectors()
                
        if parallel:
            return DoubleEndedConnectPar(min1, min2, pot, mindist, database, **kwargs)
        else:
            return DoubleEndedConnect(min1, min2, pot, mindist, database, **kwargs)

    #
    #the following functions are used only for the GUI
    #

    def draw(self, coords, index):
        """
        tell the gui how to represent your system using openGL objects
        
        Parameters
        ----------
        coords : array
        index : int
            we can have more than one molecule on the screen at one time.  index tells
            which one to draw.  They are viewed at the same time, so they should be
            visually distinct, e.g. different colors.  accepted values are 1 or 2        
        """
        raise NotImplementedError

    def smooth_path(self, images):
        """return a smoothed path between two configurations.

        used for movies
        
        See Also
        --------
        pygmin.landscape.smoothPath
        """
        raise NotImplementedError
    
    def createNEB(self):
        """ """
        raise NotImplementedError
    
    def load_coords_pymol(self, coordslist, oname, index=1):
        """load the coords into pymol
        
        the new object must be named oname so we can manipulate it later
                        
        Parameters
        ----------
        coordslist : list of arrays
        oname : str
            the new pymol object must be named oname so it can be manipulated
            later
        index : int
            we can have more than one molecule on the screen at one time.  index tells
            which one to draw.  They are viewed at the same time, so they should be
            visually distinct, e.g. different colors.  accepted values are 1 or 2
        
        Notes
        -----
        the implementation here is a bit hacky.  we create a temporary xyz file from coords
        and load the molecule in pymol from this file.  
        """
        #pymol is imported here so you can do, e.g. basinhopping without installing pymol
        import pymol 

        #create the temporary file (.xyz or .pdb, or whatever else pymol can read)
        #note: this is the part that will be really system dependent.        
        f = tempfile.NamedTemporaryFile(mode="w", suffix=".xyz")
        fname = f.name
        #write the coords into file
        for coords in coordslist:
            write_xyz(f, coords, title=oname)
        f.flush()
                
        #load the molecule from the temporary file
        pymol.cmd.load(fname)
        
        #get name of the object just create and change it to oname
        objects = pymol.cmd.get_object_list()
        objectname = objects[-1]
        pymol.cmd.set_name(objectname, oname)
        
        #here you might want to change the representation of the molecule, e.g.
        # >>> pymol.cmd.hide("everything", oname)
        # >>> pymol.cmd.show("spheres", oname)
        
        #set the color according to index
        if index == 1:
            pymol.cmd.color("red", oname)
        else:
            pymol.cmd.color("gray", oname)



if __name__ == "__main__":
    mysys = BaseSystem()
    mysys.get_potential()
    mysys.get_basinhopping()
