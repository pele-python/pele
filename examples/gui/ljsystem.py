from PyQt4 import QtGui
import NewLJ
import sys
import numpy as np
from storage import savenlowest
import time
from NEB import NEB
        
class LJSystem:
    def __init__(self):
        dlg = NewLJDialog()
        dlg.exec_()
        self.natoms = dlg.natoms()
        self.nsave = dlg.nsave()
        if dlg.result() == QtGui.QDialog.Rejected:
            raise BaseException("Aborted parameter dialog")
        self.storage = savenlowest.SaveN(self.nsave)
        
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
    
    def draw(self, coordslinear):
        from OpenGL import GL,GLUT
        coords = coordslinear.reshape(coordslinear.size/3, 3)
        com=np.mean(coords, axis=0)                  
        for xx in coords:
            x=xx-com
            GL.glPushMatrix()            
            GL.glTranslate(x[0],x[1],x[2])
            GLUT.glutSolidSphere(0.5,30,30)
            GL.glPopMatrix()
    
    def Align(self, coords1, coords2):
        from mindist.minpermdist_stochastic import minPermDistStochastic as minpermdist
        dist, X1, X2 = minpermdist( coords1, coords2, niter = 100 )
        return X1, X2
    
    def createNEB(self, coords1, coords2):
        import potentials.lj as lj
        return NEB.NEB(coords1, coords2, lj.LJ(), k = 100. ,nimages=20)

       

class NewLJDialog(QtGui.QDialog,NewLJ.Ui_DialogLJSetup):
    def __init__(self):
        QtGui.QDialog.__init__(self)
        self.setupUi(self)
    def natoms(self):
        return int(self.lineNatoms.text())
    def nsave(self):
        return int(self.lineNsave.text())
        
if __name__ == "__main__":
    import gui.run
    gui.run.run_gui(LJSystem)