import sys

from PyQt4 import QtGui

from pele.systems import LJCluster
from pele.gui import run_gui

"""
start a gui for a lennard jones cluster.

All that is really needed to start a gui is define a system and call run_gui

    system = LJCluster(13)
    run_gui(system)

"""

if __name__ == "__main__":
    # create a pop up window to get the number of atoms
    app = QtGui.QApplication(sys.argv)
    dialog = QtGui.QInputDialog()
    dialog.setLabelText("number of atoms")
    dialog.setWindowTitle("Create new Lennard-Jones system")
    dialog.setInputMode(1)
    dialog.setIntMinimum(2)
    dialog.setIntValue(13)
    dialog.exec_()
    if not dialog.result():
        sys.exit()
    natoms = dialog.intValue()

    # create the system and start the gui
    # (note: since the application is already started we need to pass it to run_gui)
    system = LJCluster(natoms)
    run_gui(system, application=app)
