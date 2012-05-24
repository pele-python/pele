from PyQt4 import QtGui
import NewLJ
import sys
import numpy as np
from storage import savenlowest
import time
from NEB import NEB
        
class BLJSystem:
    def __init__(self):
        dlg = NewLJDialog()
        dlg.exec_()
        self.natoms = dlg.natoms()
        self.ntypeA = int(0.8*self.natoms)
        self.nsave = dlg.nsave()
        if dlg.result() == QtGui.QDialog.Rejected:
            raise BaseException("Aborted parameter dialog")
        self.storage = savenlowest.SaveN(self.nsave)
        
    def createBasinHopping(self):
        import basinhopping as bh
        import potentials.ljpshiftfast as lj
        from take_step import random_displacement
        coords = np.random.random(3 * self.natoms)
        potential = lj.LJpshift(self.natoms, self.ntypeA)
        self.sigBB = potential.BB.sig
        step = random_displacement.takeStep(stepsize=0.5)
        opt = bh.BasinHopping(coords,potential,
                          temperature=1., takeStep=step)
        return opt
    
    def draw(self, coordslinear, index):
        # index = 1 or 2
        from OpenGL import GL,GLUT
        coords = coordslinear.reshape(coordslinear.size/3, 3)
        com=np.mean(coords, axis=0)                  
        size = 0.5
        if index == 1:
            color = [0.65, 0.0, 0.0, 1.]
        else:
            color = [0.00, 0.65, 0., 1.]
        GL.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE, color)
        for i,xx in enumerate(coords):
            if i == self.ntypeA: 
                size *= 0.88 #this should be dependent on lj parameters
                if index == 1:
                    color = [0.25, 0.00, 0., 1.]
                else:
                    color = [0.00, 0.25, 0., 1.]
                GL.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE, color)
            x=xx-com
            GL.glPushMatrix()            
            GL.glTranslate(x[0],x[1],x[2])
            GLUT.glutSolidSphere(size,30,30)
            GL.glPopMatrix()
    
    def Align(self, coords1, coords2):
        from mindist.minpermdist_stochastic import minPermDistStochastic as minpermdist
        self.permlist = [range(self.ntypeA), range(self.ntypeA, self.natoms)] #permutable atoms
        dist, X1, X2 = minpermdist( coords1, coords2, niter = 100, permlist = self.permlist )
        return X1, X2
    
    def createNEB(self, coords1, coords2):
        import potentials.ljpshiftfast as lj
        return NEB.NEB(coords1, coords2, lj.LJpshift(self.natoms, self.ntypeA), k = 100. ,nimages=20)

       

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
    gui.run.run_gui(BLJSystem)
