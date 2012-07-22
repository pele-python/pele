from PyQt4 import QtGui
import NewLJ
import sys
import numpy as np
from pygmin.storage import savenlowest
import time
from pygmin.NEB import NEB,dimer,tstools
import pygmin.potentials.lj as lj
                
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
        import pygmin.basinhopping as bh
        import pygmin.potentials.lj as lj
        from pygmin.takestep import displace
        coords = np.random.random(3 * self.natoms)
        potential = lj.LJ()
        step = displace.RandomDisplacement(stepsize=0.5)
        opt = bh.BasinHopping(coords,potential,
                          temperature=1., takeStep=step)
        return opt
    
    def draw(self, coordslinear, index):
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
        from pygmin.mindist.minpermdist_stochastic import minPermDistStochastic as minpermdist
        dist, X1, X2 = minpermdist( coords1, coords2, niter = 100 )
        return X1, X2
    
    def createNEB(self, coords1, coords2):
        return NEB.NEB(coords1, coords2, lj.LJ(), k = 100. ,nimages=20)
    
    def findTS(self, coords):
        pot = lj.LJ()
        tau = np.random.random(coords.shape) - 0.5
        x0,E,tau = dimer.findTS(pot, coords+5e-1*tau, tau)
        m1,m2 = tstools.minima_from_ts(pot.getEnergyGradient, x0, tau, displace=1e-2)
        print "Energies: ", m1[1],E,m2[1]
        return [x0,E],m1,m2
    
        
       

class NewLJDialog(QtGui.QDialog,NewLJ.Ui_DialogLJSetup):
    def __init__(self):
        QtGui.QDialog.__init__(self)
        self.setupUi(self)
    def natoms(self):
        return int(self.lineNatoms.text())
    def nsave(self):
        return int(self.lineNsave.text())
        
if __name__ == "__main__":
    import pygmin.gui.run as gr
    gr.run_gui(LJSystem)
