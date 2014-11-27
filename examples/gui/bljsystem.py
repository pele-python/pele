"""
start a gui for a binary lennard jones cluster.

All that is really needed to start a gui is define a system and call run_gui

    system = BLJCluster(natoms, ntypeA)
    run_gui(system)

"""
import sys

from PyQt4 import QtGui

from pele.systems import BLJCluster
from pele.gui import run_gui

from _blj_dialog import Ui_DialogLJSetup as UI


class BLJDialog(QtGui.QDialog):
    def __init__(self):
        QtGui.QDialog.__init__(self)
        self.ui = UI()
        self.ui.setupUi(self)
        self.setWindowTitle("Create binary Lennard-Jones system")
        self.natoms = None
#        self.ui.buttonBox.Ok.setDefault(True)
#        self.ui.buttonBox.Ok.setDefault(True)
    
    def get_input(self):
        self.natoms = int(self.ui.lineEdit_natoms.text())
        self.ntypeA = int(self.ui.lineEdit_ntypeA.text())
        self.sigAB = float(self.ui.lineEdit_sigAB.text())
        self.epsAB = float(self.ui.lineEdit_epsAB.text())
        self.sigBB = float(self.ui.lineEdit_sigBB.text())
        self.epsBB = float(self.ui.lineEdit_epsBB.text())
        self.sigAA = 1.
        self.epsAA = 1.

    def on_buttonBox_accepted(self):
        self.get_input()
        self.close()

    def on_buttonBox_rejected(self):
        self.close()

if __name__ == "__main__":
    # create a pop up window to get the number of atoms
    app = QtGui.QApplication(sys.argv)
    dialog = BLJDialog()
    dialog.exec_()
    
    if dialog.natoms is None:
        sys.exit()

    print dialog.ntypeA, "A atoms interacting with eps", dialog.epsAA, "sig", dialog.sigAA
    print dialog.natoms - dialog.ntypeA, "B atoms interacting with eps", dialog.epsBB, "sig", dialog.sigBB
    print "The A and B atoms interact with eps", dialog.epsAB, "sig", dialog.sigAB
    
    # create the system and start the gui
    # (note: since the application is already started we need to pass it to run_gui)
    system = BLJCluster(dialog.natoms, dialog.ntypeA,
                        sigAB=dialog.sigAB,
                        epsAB=dialog.epsAB,
                        sigBB=dialog.sigBB,
                        epsBB=dialog.epsBB,
                        )
    run_gui(system, application=app)
