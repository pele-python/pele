.. _optim2params:

The Parameter tree and OPTIM keywords
-------------------------------------

This page describes how to translate an OPTIM odata file into the pele
:ref:`Parameter <system_class>` tree.  If a keyword is not mentioned, then it is either not relevant, or
I just haven't gotten to it yet.

The documentation for OPTIM keywords can be found here, `<http://www-wales.ch.cam.ac.uk/OPTIM.doc/node4.html>`_

Generally shortened names are used for the Parameter tree, e.g. `lowestEigenvectorQuenchParams` instead
of `params.double_ended_connect.local_connect_params.tsSearchParams.lowestEigenvectorQuenchParams`.  See
:ref:`here <system_class>` for a full description of the Parameter tree.

NEWCONNECT NConMax NTriesMax ImageDensity IterDensity ImageMax ImageIncr RMStol::

  double_ended_connect.niter = NConnMax
  (local_connect.NEBattempts = NTriesMax not passable)
  NEBparams.image_density =  ImageDensity
  NEBparams.iter_density =  IterDensity
  NEBparams.max_images =  ImageMax
  (ImageIncr not implemented)
  NEBparams.NEBquenchParams.tol = RMStol

NEBK nebk::

  NEBparams.k = nebk

MAXERISE x1 x2::

  (x1 is for normal quenches)
  structural_quench_params.maxErise = x1
  lowestEigenvectorQuenchParams.maxErise = x2

BFGSTS nevs ncgmax1 ncgmax2 ceig nevl::

  lowestEigenvectorQuenchParams.nsteps = nevs
  tsSearchParams.nsteps_tangent1 = ncgmax1
  tsSearchParams.nsteps_tangent2 = ncgmax2
  lowestEigenvectorQuenchParams.tol = ceig
  (nevl not used)

BFGSMIN gmax::

  structural_quench_params.tol = gmax

MAXSTEP x1 x2::

  (x1 is for eigenvector-following, but is it used for other things?)
  (x2 ? not documented)

MAXMAX x::

  tsSearchParams.max_uphill_step = x

BFGSCONV gmax::

  tsSearchParams.tangentSpaceQuenchParams.tol = gmax (I think)

PUSHOFF x::

  local_connect_params.pushoff_params.stepmin = x

STEPS n::

  ?

BFGSSTEPS n::

  (number of steps for standard quenches)
  structural_quench_params.nsteps = n

MAXBFGS x1 x2 x3 x4::

  (max step length)
  (x1 for normal quenches)
  structural_quench_params.maxstep = x1 (maybe other places also?)
  lowestEigenvectorQuenchParams.maxstep = x2
  (x3 for putting structures in closest coincidence with ministd)
  NEBparams.NEBquenchParams.maxstep = x4

UPDATES mupdate1 mupdate2 mupdate3 mupdate4::

  (M for lbfgs)
  (mupdate1 is for standard quenches)
  structural_quench_params.M = mupdate1 (maybe other places also?)
  (??? mupdate2 is currently ignored by mind)
  lowestEigenvectorQuenchParams.M = mupdate3
  NEBparams.NEBquenchParams.M = mupdate4

CHECKCHIRALITY::
not implemented, but should be


RANROT n::

  (niter in the stochastic minpermdist routine.  not in parameter tree yet)

The Parameter tree and GMIN keywords
-------------------------------------
EDIFF econv::

  params.database.accuracy = econv

SLOPPYCONV::
  no equivalent, see TIGHTCONV.  Also, see :meth:`.BaseSystem.get_compare_exact`

TIGHTCONV cgmax::

  params.structural_quench_params.tol = cgmax

UPDATES nup::

  params.structural_quench_params.M = nup

MAXBFGS max::

  params.structural_quench_params.maxstep = nup

