"""
.. currentmodu:: pele.gui

pele GUI (`pele.gui`)
--------------------------
This module contains all the necessary components for running the gui.  

.. autosummary::
    :toctree: generated/

    run_gui

Simply initialize your system and pass it to run_gui::

    from pele.systems import LJCluster
    from pele.gui import run_gui
    mysystem = LJCluster(13)
    run_gui(mysystem)

if you pass a database file name it will connect to an existing database
or create a new one at that location::

    run_gui(mysystem, db="database.sqlite")


Preparing my System for use in the gui
++++++++++++++++++++++++++++++++++++++
If you have written your own system class and want to run it in the gui, there are
a few additional member functions must be defined::

    mysystem.draw(coords, index)
    mysystem.smooth_path(images)

Both of these are used for displaying your system in a 3D OpenGL renderer.
See :ref:`BaseSystem <system_class>` and existing derived classes like LJCluster
and BLJCluster for more information and examples of how to implement these.

"""
from __future__ import absolute_import






from .run import run_gui
#from ui.mplwidget import MPLWidget

