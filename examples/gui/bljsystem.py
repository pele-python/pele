from PyQt4 import QtGui
import NewLJ
import sys
import numpy as np
#from storage import savenlowest
import time
#from NEB import NEB

from pygmin import gui
import pygmin.basinhopping as bh
import pygmin.potentials.ljpshiftfast as lj
from pygmin.takestep import displace
from pygmin.mindist import ExactMatchCluster, MinDistWrapper
from pygmin.mindist import minPermDistRanRot as minpermdist
from ljsystem import LJSystem


class CompareExact(object):
    def __init__(self, accuracy = 1e-2, **kwargs):
        self.kwargs = kwargs
        self.compare_exact = ExactMatchCluster(**kwargs)
    
    def __call__(self, m1, m2):
        ret =  self.compare_exact(m1.coords, m2.coords)
        if( not ret):
            print "2 minima within enery accuracy have different coordinates", m1.energy, m2.energy
        return ret


class BLJSystem(LJSystem):
    def __init__(self):
#        dlg = NewLJDialog()
#        dlg.exec_()
#        self.natoms = dlg.natoms()
#        self.nsave = dlg.nsave()
#        if dlg.result() == QtGui.QDialog.Rejected:
#            raise BaseException("Aborted parameter dialog")
        
        super(BLJSystem, self).__init__()
        print self.natoms
        self.ntypeA = int(0.8*self.natoms)

    def set_database(self, database):
        permlist = [range(self.ntypeA), range(self.ntypeA, self.natoms)]
        compare = CompareExact(permlist=permlist)
        database.compareMinima = compare
        gui.GUISystem.set_database(self, database)


    def create_potential(self):
        potential = lj.LJpshift(self.natoms, self.ntypeA)
        self.sigBB = potential.BB.sig
        return potential
    
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

    def create_mindist_object(self): 
        self.permlist = [range(self.ntypeA), range(self.ntypeA, self.natoms)] #permutable atoms
        return MinDistWrapper(minpermdist, permlist=self.permlist)

       
from ljsystem import NewLJDialog
#class NewLJDialog(QtGui.QDialog,NewLJ.Ui_DialogLJSetup):
#    def __init__(self):
#        QtGui.QDialog.__init__(self)
#        self.setupUi(self)
#    def natoms(self):
#        return int(self.lineNatoms.text())
#    def nsave(self):
#        return int(self.lineNsave.text())
        
if __name__ == "__main__":
    import pygmin.gui.run as gr
    gr.run_gui(BLJSystem)
