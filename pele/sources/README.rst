Pele's C++ Backend
------------------

Much of the backend of pele is written in C++ in order for the calculations to
be as fast as possible.  Primarily this means the potentials, some of the
optimizers, and a Monte Carlo and Parallel temperring routine (coming soon).  
Initially we had done this with previous implementations in various languages
(c, fortran, cython).  These were tied together in Python, which had the drawback that 
a fortran minimizer needed to go through python to get access to a Fortran potential function.
The current C++ unified framework allows us to avoid the slow python gatekeeper.
In this framework a C++ potential and optimizer are created in python.  The C++
potential is passed directly to the C++ optimizer and the entire minimization
is done without need to reference any Python functions.
All of the interaction between Python and C++ is done through Cython.

pele::Array 
----------- 

The pele::Array duplicates the functionality of an std::vector or a numpy
array.  The main reason for using Array is to have a pure C++ array class that
can wrap a data array (c array) without owning it.  For example, the data in a numpy array
is wrapped by a pele::Array in cython and passed to the C++ function.  This avoids
needless copying of potentially large data arrays.

pele::Arrays are used everywhere that arrays need to be passed around,
especially for functions that are exposed to cython.  For all actions with an
Array (copy constructor, assignment operator, etc) the default behavior
follows python conventions and wraps the data rather than copy it.  A
pele::Array uses reference counting (similar to std::shared_ptr) to ensure that
data in an Array is freed only when it is no longer needed by any Arrays.

BasePotential
-------------

The C++ class BasePotential is the base class for all potential objects.  Every
c++ potential must derive from BasePotential and, at a minimum, overload the
function get_energy().  It is recommended to overload get_energy_gradient() as
well.  The energy functions can be computed any way you want, but for simple
potentials life is easier if you separate the interaction from the looping over
atoms.  You write your own interaction and use an existing looping class like
pele::simple_pairwise_potential.  With C++ templates and inline functions you
won't sacrafice any speed by structuring it like this.  See the Lennard-Jones 
potential `lj.h` for an example of how to set this up.

Eclipse CDT
-----------
Eclipse is a bit buggy and you have to do a bit of fiddling to get it to
recognize c++ and c++11 features.  I've found this answer to be the simplest
way of getting it working.  http://stackoverflow.com/a/13549029/3307093
