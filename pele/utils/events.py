"""
adapted from http://code.activestate.com/recipes/577980-improved-signalsslots-implementation-in-python/

A signal/slot implementation

File:    signal.py
Author:  Thiago Marcos P. Santos
Author:  Christopher S. Case
Author:  David H. Bronke
Created: August 28, 2008
Updated: December 12, 2011
License: MIT

"""

import inspect
from weakref import WeakSet, WeakKeyDictionary

__all__ = ["Signal"]


class Signal(object):
    """ class for signal slot concept

    Example
    -------

    A simple example for a callback is
    >>> event = Signal()
    >>> event.connect(mfunc)
    >>> # raise the signal
    >>> event("hello")
    >>>
    >>> # functions can be disconnected
    >>> even.disconnect(myfunc)

    Since weak references are used, care has to be taken with object functions

    >>> obj = MyClass()
    >>> event.connect(obj.myfunc) # works
    >>> event.connect(MyClass().myfunc) # will not work

    The second example for member functions will not work since the Signal class
    uses weakref and therefore does not increase the reference counter. MyClass()
    only exists for the time of the function call and will be deleted afterwards
    and the weakref will become invalid.

    """
    
    def __init__(self):
        self._functions = WeakSet()
        self._methods = WeakKeyDictionary()

    def __call__(self, *args, **kargs):
        """ raise the event """
        # Call handler functions
        for func in self._functions:
            func(*args, **kargs)

        # Call handler methods
        for obj, funcs in self._methods.items():
            for func in funcs:
                func(obj, *args, **kargs)

    def connect(self, slot):
        """ connect a function / member function to the signal """
        if inspect.ismethod(slot):
            if slot.__self__ not in self._methods:
                self._methods[slot.__self__] = set()

            self._methods[slot.__self__].add(slot.__func__)

        else:
            self._functions.add(slot)

    def disconnect(self, slot):
        """ disconnect a function from the signal """
        if inspect.ismethod(slot):
            if slot.__self__ in self._methods:
                self._methods[slot.__self__].remove(slot.__func__)
        else:
            if slot in self._functions:
                self._functions.remove(slot)

    def clear(self):
        """ remove all callbacks from the signal """
        self._functions.clear()
        self._methods.clear()

