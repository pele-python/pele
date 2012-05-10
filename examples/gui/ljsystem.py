from PyQt4 import QtGui
import NewLJ
import sys

class LJSystem:
    def __init__(self):
        dlg = NewLJDialog()
        dlg.exec_()
        self.natoms = dlg.natoms()
        self.nsave = dlg.nsave()
        if dlg.result() == QtGui.QDialog.Rejected:
            raise BaseException("Aborted parameter dialog")
        
    def createBasinHopping(self):
        import numpy as np
        import basinhopping as bh
        import potentials.lj as lj
        from take_step import random_displacement
        coords = np.random.random(3 * self.natoms)
        potential = lj.LJ()
        step = random_displacement.takeStep(stepsize=0.5)
        opt = bh.BasinHopping(coords,potential,
                          temperature=1., takeStep=step.takeStep)
        return opt
       

class NewLJDialog(QtGui.QDialog,NewLJ.Ui_DialogLJSetup):
    def __init__(self):
        QtGui.QDialog.__init__(self)
        self.setupUi(self)
    def natoms(self):
        return int(self.lineNatoms.text())
    def nsave(self):
        return int(self.lineNsave.text())
        
if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    a = LJSystem()
    print a.natoms
    #print a.nsave