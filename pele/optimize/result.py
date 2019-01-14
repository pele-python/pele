from __future__ import print_function
# Result object copied from scipy 0.11
__all__ = ['Result']


class Result(dict):
    """A container for the return values of an optimizer

    Attributes
    ----------
    coords : ndarray
        The solution of the optimization.
    success : bool
        Whether or not the optimizer exited successfully.
    status : int
        Termination status of the optimizer. Its value depends on the
        underlying solver. Refer to `message` for details.
    message : str
        Description of the cause of the termination.
    energy : ndarray
        energy at the solution
    grad : ndarray
        gradient at the solution
    nfev : int
        Number of evaluations of the function or gradient
    nit : int
        Number of iterations performed by the optimizer.

    Notes
    -----
    There may be additional attributes not listed above depending of the
    specific solver. 
    
    Also, since this class is essentially a subclass of dict
    with attribute accessors, one can see which attributes are available
    using the `keys()` method.
    """

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name)

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


    def __repr__(self):
        if self.keys():
            m = max(list(map(len, self.keys()))) + 1
            return '\n'.join([k.rjust(m) + ': ' + repr(v)
                              for k, v in self.items()])
        else:
            return self.__class__.__name__ + "()"

    def __getitem__(self, i):
        """
        April 26, 2013
             
        this overloaded function exists only for compatibility with the old quenchers.
        It should be removed at some point in the future
        """
        if isinstance(i, slice):
            mylength = 5
            vals = tuple([self[j] for j in range(*i.indices(mylength))])
            return vals

        if i in range(4):
            maplist = {0: "coords",
                       1: "energy",
                       2: "rms",
                       3: "nfev",
            }
            i = maplist[i]
        elif i == 4:
            return self
        return dict.__getitem__(self, i)

    def __iter__(self):
        """
        April 26, 2013
             
        this overloaded function exists only for compatibility with the old quenchers.
        It should be removed at some point in the future
        """
        return iter((self.coords, self.energy, self.rms, self.nfev, self))


if __name__ == "__main__":
    import numpy as np

    res = Result()
    res.coords = np.array([0])
    res.energy = 1.
    res.rms = 1e-4
    res.nfev = 100
    print(dir(res))
    print(res[0])
    print(res[1])
    print(res[4])
    print("slice", res[:2])
    x, e = res[:2]
    print("unpack slice", x, e)
    a, b, c, d, f = res
    print("done")
    print(a, b, c, d)
