PyQt5 GUI audit
==============

Expected behavior
-----------------

* Launch the Lennard-Jones example with a new database or open an existing one;
  render structures with working rotation, zoom, and window resizing.
* Switch databases without retaining old rows, selections, graph windows, or
  background jobs. A failed open must preserve the current database and views.
* Display minima and transition states with numeric energy/ID sorting and the
  configured list limit; reflect database additions, removals, and merges.
* Select a minimum in the main list, two independent connection endpoints,
  or a transition state and its adjacent minima. Graph picks must work even
  when a minimum is outside the displayed list limit.
* Keep selection-dependent actions unavailable until their inputs exist.
  Right-clicking empty lists and stopping an idle calculation must be harmless.
* Export coordinates, preserving their values, and cancel file dialogs safely.
* Start basin-hopping with a valid step count, update the worker count, retain
  all returned minima, stop responsively, and clean up on exit.
* Align structures without changing database coordinates; connect/reconnect
  selected minima, display completed paths, and report worker failures.
* Run Connect All, stop/restart, toggle result views before and after completion,
  and close windows without leaving background processes running.
* Initialize, optimize, reset, save, and reload a nudged elastic band; inspect
  its energies and frames; select and refine transition-state candidates.
* Browse normal modes for minima and transition states, show mode animations,
  and handle an empty selection or invalid amplitude without a traceback.
* Generate, perturb, quench, and save structures in the takestep explorer.
* Edit nested numeric, string, and Boolean parameters with changes reaching
  the system; reject invalid numeric edits without corrupting the value.
* Draw empty and populated connectivity/disconnectivity graphs; pick minima,
  zoom, redraw, adjust levels, highlight minimum-energy paths, and export plots.
* Color graph nodes by committor and mean first passage time using all relevant
  transition states and the selected target, with consistent repeated results.
* Compute thermodynamic information asynchronously; compute heat capacity and
  rates with valid inputs; report empty-data/invalid-input errors in the dialog;
  cancel workers when dialogs close.
* Close all analysis windows and exit cleanly; reopen tools against current data.
* Open About. Optional PyMOL and external OPTIM integration should work when
  those programs and a compatible system configuration are available.

Audit outcome
-------------

The audit added regression checks for failures in database switching, deletion,
merging, limited lists, empty selections, coordinate export, parameter editing,
NEB loading/frame selection, normal modes, animation, graph calculations,
worker result transfer, cancellation, restart, and input validation.

Desktop verification used LJ13: basin-hopping returned five minima; the actual
OpenGL framebuffer rendered the cluster; normal modes and heat capacity opened
and calculated successfully; analysis windows and workers closed cleanly.
Visual inspection also caught and corrected excess blank space above the
heat-capacity plot. Graph interaction and plot export were checked separately.

Scientific workflow checks used LJ7 for a real basin-hopping/connect sequence
(four minima found, successful connection with two transition states) and LJ3
for interactive NEB, normal-mode animation/energies, and displacement/quenching.

The optional PyMOL and external OPTIM integrations were not exercised. Parameter
context-menu Add/Delete actions were unfinished placeholders; they are now
disabled rather than appearing to perform an operation. Existing values remain
editable. Forced shutdown may discard work that a stuck process has not yet
transferred; normal returned results are drained before cleanup.

Verification
------------

The regression suite uses real PyQt5 widgets, SQLite databases, scientific
calculations, and small multiprocessing jobs. Headless tests replace GLUT
initialization only where a desktop display is unavailable::

  MPLCONFIGDIR=/tmp/pele-gui-matplotlib QT_QPA_PLATFORM=offscreen \
    OMP_NUM_THREADS=1 python -m pytest pele/gui/tests -q

Core checks, including the shared thermodynamic and rate calculations::

  OMP_NUM_THREADS=1 python -m pytest pele/ -q

Final run in the ``pele-meson`` Python 3.13 environment: 661 tests and 19
subtests passed. The run emitted 220 existing deprecation/syntax warnings.

The checklist describes expected behavior; passing regression tests does not
certify every potential/system implementation or optional external integration.
