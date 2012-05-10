from PyQt4 import QtGui
import NewLJ
import sys
import numpy as np
        
class LJSystem:
    def __init__(self):
        dlg = NewLJDialog()
        dlg.exec_()
        self.natoms = dlg.natoms()
        self.nsave = dlg.nsave()
        if dlg.result() == QtGui.QDialog.Rejected:
            raise BaseException("Aborted parameter dialog")
        
    def createBasinHopping(self):
        import basinhopping as bh
        import potentials.lj as lj
        from take_step import random_displacement
        coords = np.random.random(3 * self.natoms)
        potential = lj.LJ()
        step = random_displacement.takeStep(stepsize=0.5)
        opt = bh.BasinHopping(coords,potential,
                          temperature=1., takeStep=step.takeStep)
        return opt
    
    def draw(self, coords):
        from OpenGL import GL,GLUT
        com=np.mean(coords, axis=0)                  
        for xx in coords:
            x=xx-com
            GL.glPushMatrix()            
            GL.glTranslate(x[0],x[1],x[2])
            GLUT.glutSolidSphere(0.5,30,30)
            GL.glPopMatrix()
       
       

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