.. _tutorial_system_class:

Using the system class
----------------------

All of the tools in pele are modular, and can be used without a system class,
but seting up a system class can make your life a *lot* easier.  If you set it up
once then access to all the tools in pele is simply one function call away.  This
tutorial is designed to show how use the system class to access the most common
features of pele.

.. note::
  This tutorial mirrors an example in the examples folder.

This tutorial uses the LJCluster system class as an example.  We
start by initializing the system class.  Different system classes
may take different parameters for the constructor.  LJCluster simply
takes the number of atoms.  We'll look at clusters of 13 atoms::

  from pele.systems import LJCluster
  natoms = 13
  system = LJCluster(natoms)

We start by generating a random configuration::

  coords = system.get_random_configuration()

We want to know the energy of that configuration, so we use the system
class to create a :ref:`potential object <potentials_description>`.  This
potential object has functions for computing the energy, gradient, and
Hessian::

  potential = system.get_potential()
  energy = potential.getEnergy(coords)
  print "the energy of the random configuration is", energy

Next we minimize that configuration to it's nearest local minimum. To do 
that we use the system class to generate an object to do the :ref:`minimization <optimize_module>`
for us::

  quencher = system.get_minimizer()
  ret = quencher(coords)
  newcoords = ret.coords
  newenergy = ret.energy
  print "after quenching, the energy is", newenergy

Next we will create a :ref:`database <database>` to store the newly
generated minimum::

  db = system.create_database()
  minimum1 = db.addMinimum(newenergy, newcoords)

.. note::
  This creates a database in memory.  If you want to save the results to a file
  you would use::

    db = system.create_database("lj.sqlite")

Let's generate a second random minimum and add it to the database.  Now we 
use a shortcut from the system class::

  ret = system.get_random_minimized_configuration()
  print "a second minimum has energy", ret[1]
  e2, coords2 = ret[1], ret[0]
  minimum2 = db.addMinimum(e2, coords2)

Now we'll do a short basinhopping run to find the global minimum and build up
the database of minima::

  bh = system.get_basinhopping(database=db)
  niter = 20
  bh.run(niter)
  print "the lowest energy found after", niter, " basinhopping steps is", db.minima()[0].energy


Print the energies of all the minima we've found::

  print "the minima in the database have energies"
  for minimum in db.minima():
      print "  ", minimum.energy

Next, lets find the :ref:`minimum distance <structure_alignment>` (a.k.a. mindist) between the two lowest
minima::

  m1, m2 = db.minima()[:2]
  mindist = system.get_mindist()
  dist, coords1, coords2 = mindist(m1.coords, m2.coords)
  print "the minimum distance between the two lowest minima is", dist

Now lets do a :ref:`double ended connect <landscape_description>` run.  This
finds a connected series of minima and transition states between two minima::

  connect = system.get_double_ended_connect(m1, m2, db)
  connect.connect()
  mints, S, energies = connect.returnPath()
  nts = (len(mints) - 1)/2
  print "found a connection with", nts, "transition states"

Finally, lets connect all of the minima in the database to the lowest minimum::

  print "now connecting all the minima to the lowest energy minimum"
  from pele.landscape import ConnectManager
  manager = ConnectManager(db, strategy="gmin")
  for i in xrange(db.number_of_minima()-1):
      print "connecting minima with id's", m1._id, m2._id
      m1, m2 = manager.get_connect_job()
      connect = system.get_double_ended_connect(m1, m2, db)
      connect.connect()


And we'll end by printing out some information about what is in the database::

  print "database summary:"
  print "    ", len(db.minima()), "minima"
  print "    ", len(db.transition_states()), "transition states"
