.. _global_optimization:

Global Optimization
===================
 
The main algorithm for global optimization is basin hopping. The core object for a
basin hopping run is a BasinHopping object. The step taking and various other customization
routines (e.g. storage, acceptance criterion, ...) can be attached to this object to customize
the behaviour of the basin hopping procedure.

.. currentmodule:: pele.basinhopping

.. autosummary::
    :toctree: generated/

    BasinHopping

::

  import numpy as np
  import pele.potentials.lj as lj
  import pele.basinhopping as bh
  from pele.takestep import displace

  natoms = 12
 
  # random initial coordinates
  coords=np.random.random(3*natoms)
  potential = lj.LJ()

  step = displace.RandomDisplacement(stepsize=0.5)

  opt = bh.BasinHopping(coords, potential, takeStep=step)
  opt.run(100)


.. automodule:: pele.takestep





