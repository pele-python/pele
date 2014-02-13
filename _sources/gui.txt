.. _gui:

pele GUI
========
pele is primarily designed to be used as a library, but there is an extensive
GUI as well.  The gui can be used to explore almost everything implemented in
pele.  An example workflow might be:

1. Run basinhopping to find the global minimum structure.

2. See a list of minima found and view the structure in 3D using openGL.

.. image:: gui_main.png
  :height: 300

3. View the normal modes of each structure.

.. image:: gui_normalmodes.png
  :height: 300

4. Find transition states between minima using the double ended connect algorithm

.. image:: gui_connect_all.png
  :height: 300

5. Compare the structures of minima.  

.. image:: gui_ts.png
  :height: 300

6. See a list of all known transition
   states, view the 3d structure, and look at the normal modes.

7. Create an interactive disconnectivity graph.  Clicking on a minimum in the
   disconnectivity graph will bring up the minimum in the main window for 3D
   viewing and analysis.

.. image:: gui_dgraph.png
  :height: 300

8. Run the nudged elastic band algorithm interactively.

.. image:: gui_NEB.png
  :height: 300


.. automodule:: pele.gui

