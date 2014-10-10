"""
Global optimization
-------------------

We primarily use the global optimization method :class:`.BasinHopping` because it has been
shown to be very effective on the very high dimensional, smooth landscapes that we work with.

BasinHopping is iterative with each cycle composed of the following
features

1) random perturbation of the coordinates

2) local minimization

3) accept or reject the new coordinates based on the minimized function
   value


Let's do an example where we find run :class:`.BasinHopping` on a Lennard-Jones cluster (:class:`.LJ`).

.. note::
    This example involves quite a number of steps.  Normally these are done automatically
    in functions defined in the system class (:mod:`.systems`).  With the system class the
    following can be done in just a few steps::
        from pele.systems import LJCluster
        natoms = 17
        system = LJCluster(natoms)
        database = system.create_database('lj17.sqlite')
        bh = system.get_basinhopping(database)
        bh.run(10)
    The following example shows how to do this manually


We start by defining the potential and choosing a random set of starting coordinates
::
    import numpy as np
    from pele.potentials import LJ
    natoms = 17
    potential = LJ()
    x0 = np.random.uniform(-1, 1, 3*natoms)

We also need to set up the take-step class (:mod:`.takestep`).  We will use simple random
displacements of the coordinates (:class:`.RandomDisplacement`) and wrap it with the class 
:class:`.AdaptiveStepsizeTemperature` that
adjusts both the stepsize and temperature to achieve optimum results.
::
    from pele.takestep import RandomDisplacement, AdaptiveStepsizeTemperature
    displace = RandomDisplacement()
    adaptive_displacement = AdaptiveStepsizeTemperature(displace)

The last step is to set up an sqlite database (:mod:`.storage`) which we use to store the 
minima that we find during the basinhopping run.
::
    from pele.storage import Database
    database = Database("lj17.sqlite")

There is now an sqlite database named "lj17.sqlite" in the folder where this script was run.
We use `sqlalchemy` to communicate with the sqlite database from python.
The minima in the database can be accessed simply by::

    for m in database.minima():
        print m.energy

You can access other attributes of the :class:`.Minimum` as `minimum.coords` in the same way.


"""
def bh_with_system_class():
    from pele.systems import LJCluster
    natoms = 17
    system = LJCluster(natoms)
    database = system.create_database('lj17.sqlite')
    bh = system.get_basinhopping(database)
    bh.run(10)

def bh_no_system_class():
    import numpy as np
    from pele.potentials import LJ
    natoms = 17
    potential = LJ()
    x0 = np.random.uniform(-1, 1, 3*natoms)
    
    from pele.takestep import RandomDisplacement, AdaptiveStepsizeTemperature
    displace = RandomDisplacement()
    adaptive_displacement = AdaptiveStepsizeTemperature(displace)
    
    from pele.storage import Database
    database = Database("lj17.sqlite")
    
    from pele.basinhopping import BasinHopping
    bh = BasinHopping(x0, potential, adaptive_displacement, storage=database.minimum_adder)
    bh.run(10)
    
    for m in database.minima():
        print m.energy

if __name__ == "__main__":
    bh_no_system_class()