from pygmin.landscape import DoubleEndedConnect
from pygmin import basinhopping
from pygmin.storage import Database

class NotImplemented(BaseException):
    """
    The exception to return if there is a feature
    in the System is not implemented yet
    """


class BaseSystem(object):
    """
    this class defines a base class for a System object
    
    
    Description
    -----------
    The following functions need to be overloaded for running
    various routine

    
    Global Optimization
    -------------------
        get_potential : required
        get_takestep : requred
        get_random_configuration : optional
        get_compare_exact : optional
    
    Connecting Minima and Transition State Searches
    -----------------------------------------------
        get_potential : required
        get_mindist : required
        get_orthogonalize_to_zero_eigenvectors : required
        get_random_configuration : optional
        get_compare_exact : optional
    
    additionally, it's a very good idea to specify the accuracy in the database
    
    """
    def get_potential(self):
        """return the potential object"""
        raise NotImplemented

    def get_random_configuration(self):
        """a starting point for basinhopping, etc."""
        raise NotImplemented
    
    def get_compare_exact(self):
        """object that returns True if two structures are exact.
        
            true_false = compare_exact(min1, min2)
        """
        raise NotImplemented

    def get_compare_minima(self):
        """a wrapper for compare exact so in input can be in 
        Minimum Form"""
        compare = self.get_compare_exact()
        return lambda m1, m2: compare(m1.coords, m2.coords)
    
    def create_database(self, *args, **kwargs):
        """return a new database object"""
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
                compare_minima = None
                try:
                    compare_minima = self.get_compare_minima()
                    kwargs["compareMinima"] = compare_minima
                except NotImplemented:
                    pass
        except NotImplemented:
            #compareMinima is optional
            pass

        return Database(**kwargs)
    
    def get_takestep(self):
        """return the takestep object for use in basinhopping, etc."""
        raise NotImplemented

    def get_basinhopping(self, database=None, takestep=None, coords=None, **kwargs):
        """return the basinhopping object with takestep
        and accept step already implemented"""
        pot = self.get_potential()
        if coords is None:
            coords = self.get_random_configuration()
        if takestep is None:
            takestep = self.get_takestep()
        if database is None:
            database = self.load_database()
        minimum_adder = database.minimum_adder()
        bh = basinhopping.BasinHopping(coords, pot, takestep, minimum_adder, **kwargs)
        return bh

    def get_mindist(self):
        """return a mindist object that is callable with the form
        
        dist, X1new, X2new = mindist(X1, X2)
        
        Notes
        -----
        the mindist object returns returns the best alignment between two
        configurations, taking into account all global symmetries 
        """
        raise NotImplemented
    
    def get_orthogonalize_to_zero_eigenvectors(self):
        """return a which makes a vector orthogonal to the known zero
        eigenvectors (the eigenvectors with zero eigenvalues.  It should
        be callable with the form
        
            vec = orthogVec(vec, coords)
        
        """
        raise NotImplemented
    
    def get_double_ended_connect(self, min1, min2, database, **kwargs):
        """return a DoubleEndedConnect object"""
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
                lcp = kwargs["local_connect_params"] = dict()
            
            if "tsSearchParams" in lcp:
                tssp = lcp["tsSearchParams"]
            else:
                tssp = lcp["tsSearchParams"] = dict()
            
            if not "orthogZeroEigs" in tssp:
                tssp["orthogZeroEigs"] = self.get_orthogonalize_to_zero_eigenvectors()
                
        return DoubleEndedConnect(min1, min2, pot, mindist, database, **kwargs)







if __name__ == "__main__":
    mysys = BaseSystem()
    mysys.get_potential()
    mysys.get_basinhopping()