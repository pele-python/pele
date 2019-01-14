"""
.. currentmodule:: pele.storage

Database storage (`pele.storage`)
===================================

This module contains all the objects related to the Database class.
The database is the preferred method for storing minima and transition
states.  It is used in basinhopping to accumulate and keep track of
new minima as they are found.  And it is used heavily in DoubleEndedConnect
to keep track of minima and add new transition states.  The database
ensures that all minima and transition state objects are unique

.. autosummary::
   :toctree: generated/

    Database
    Minimum
    TransitionState



Minimum
-------
.. autosummary::
   :toctree: generated/

    Minimum

The Minimum object is the class uses to store a minima in the database.  Here by minimum we 
refer to a local minimum in the potential energy landscape, but it is really defined by
the energy and the coordinates.  It is also used heavily outside of the database as a way 
to unambiguously represent a minimum. Minima will compare equal only if they are the same minima.
Minima can also be hashed and used as keys or values in a dictionary::

    >>> min1 = database.addMinimum(energy1, coords1)
    >>> min2 = database.addMinimum(energy2, coords2)
    >>> min1 == min2 #  a valid comparison
    >>> neighbors = { min1:min2 } #  also fine

The above comparison relies on the fact that the database stores only unique minima (see 
below section for how it determines uniqueness).  
Thus Minimum objects should never be created manually,
you should always let the database do it.

    >>> min1 = Minimum(energy1, coords1) #  DON'T DO THIS
    >>> min2 = Minimum(energy2, coords2) #  DON'T DO THIS
    >>> min1 == min2   #  NOT OK

Minima can be access from the database in a number of ways, but most often by

    >>> database.minim() #  an iterator over minima sorted by energy

.. note::

    When a minimum is added to the database
    
        >>> minimum1 = database.addMinimum(energy, coords)
    
    the database does several checks to ensure that it is unique. First it checks
    to make sure the energy is not within `database.accuracy` of any other minimum.
    If the energy does compare similar, the database then compares the structures using
    `database.compareMinima()`.  By default this structural comparison does nothing, but
    it is a good idea for the user to set a comparison function approriate for their system.
    See :ref:`structural alignment <structure_alignment>` for more information about how
    to define this function. Only after a minimum passes these tests is it added to the 
    database as a unique Minimum.  If compares exact to an existing Minimum then that minimum
    is returned by `database.addMinimum()`. 

TransitionState
---------------
.. autosummary::
   :toctree: generated/

    TransitionState

The TransitionState object is used to store transition states in the database.
As with the Minimum class, TransitionState objects should never be created manually,
but always through the database::

    >>> ts = database.add_transtion_state(energy, coords, minimum1, minimum2, 
                                          eigenval=eigelval, eigenvec=eigenvec)

Transition states in the database are guaranteed to be unique by comparing both 
attached minima and by comparing the transition state energy.  Transition states
can be accessed from the database by::

    >>> database.transition_states #  an iterator over transition states
    >>> ts = database.getTranstionState(minimum1, minimum2)


Database 
--------

.. autosummary::
   :toctree: generated/

    Database

The Database class is what deals with creating new Minimum and TransitionState objects
and ensuring that they are unique.  It stores the objects in an SQL database on the hard drive
(or optionally in memory).
The objects are persistent in the database and exist as soon as the Database class is 
connected to the database. If any value in the objects is changed, the changes are 
automatically persistent in the database.  Database uses the pypi package SQLAlchemy to handle
the interactions with the SQL database.

To connect to a database::

    >>> db = Database() #  create a database in memory
    >>> db = Database(db="mydatabase.sqlite") 

Database will load data from "mydatabase.sqlite" if it already exists, or create a new 
database with that name if it doesn't.

.. note::

    basinhopping doesn't accept a database as a parameter, instead you should pass
    a call back function function which will add the minima as they are found.
    This function can be conveniently accessed by::
    
        >>> minima_adder = database.minimum_adder()
        >>> bh = BasinHopping(coords, potential, takestep, storage=minima_adder)

"""
from __future__ import absolute_import




from .database import *

