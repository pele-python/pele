from PyQt4 import QtGui
import NewLJ
import numpy as np
from pygmin.transition_states import NEB, dimer, tstools, InterpolatedPath
import pygmin.potentials.lj as lj
from pygmin import gui
from pygmin.mindist import ExactMatchCluster, MinDistWrapper
from pygmin.mindist import minPermDistRanRot as minpermdist
from pygmin.landscape import LocalConnect, smoothPath, DoubleEndedConnect
import pygmin.basinhopping as bh
from pygmin.takestep import displace


def compare(x1, x2):
    ret =  ExactMatchCluster(accuracy = 1e-2)(x1.coords, x2.coords)
    if( not ret):
        print "2 minima within enery accuracy have different coordinates", x1.energy, x2.energy
    return ret
                
class LJSystem(gui.GUISystem):
    def __init__(self):
        dlg = NewLJDialog()
        dlg.exec_()
        self.natoms = dlg.natoms()
        if dlg.result() == QtGui.QDialog.Rejected:
            raise BaseException("Aborted parameter dialog")
        
    def set_database(self, database):
        database.compareMinima = compare
        gui.GUISystem.set_database(self, database)
    
    def create_potential(self):
        return lj.LJ()
    
    def create_basinhopping(self):
        coords = np.random.random(3 * self.natoms)
        potential = self.create_potential()
        step = displace.RandomDisplacement(stepsize=0.5)
        opt = bh.BasinHopping(coords, potential,
                          temperature=1., takeStep=step,
                          outstream = None)
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
    
    def create_mindist_object(self):
        return MinDistWrapper(minpermdist, permlist=[range(self.natoms)], niter=10)
    
    def Align(self, coords1, coords2):
        mindist = self.create_mindist_object()
        dist, X1, X2 = mindist(coords1, coords2)
        print "after alignment the distance is", dist
        return X1, X2
    
    def createNEB(self, coords1, coords2):
        pot = self.create_potential()
        return NEB(InterpolatedPath(coords1, coords2, 30), pot, k = 100.)
    
    def zeroEigenVecs(self, coords):
        # translational eigenvectors
        x1 = np.zeros(coords.shape)
        x2 = x1.copy()
        x3 = x1.copy()
        x1.reshape(coords.size/3,3)[:,0] = 1.
        x2.reshape(coords.size/3,3)[:,1] = 1.
        x3.reshape(coords.size/3,3)[:,2] = 1.
        Rx = np.array([[ 1.,  0.,  0.],
                       [ 0.,  0.,  1.],
                       [ 0., -1.,  0.]])
        Ry = np.array([[ 0.,  0.,  1.],
                       [ 0.,  1.,  0.],
                       [-1.,  0.,  0.]])
        Rz = np.array([[ 0.,  1., 0.],
                       [-1.,  0., 0.],
                       [ 0.,  0., 1.]])
        x = coords.reshape(coords.size/3,3)
        r1 = np.dot(Rx,x.transpose()).transpose().reshape(coords.shape)
        r2 = np.dot(Ry,x.transpose()).transpose().reshape(coords.shape)
        r3 = np.dot(Rz,x.transpose()).transpose().reshape(coords.shape)

        return [x1/np.linalg.norm(x1), x2/np.linalg.norm(x2), x3/np.linalg.norm(x3), 
                r1/np.linalg.norm(r1), r2/np.linalg.norm(r2), r3/np.linalg.norm(r3)]
    
    def findTS(self, coords):
        pot = lj.LJ()
                
        ret = dimer.findTransitionState(coords+np.random.random(coords.shape)*0.01, pot, zeroEigenVecs=self.zeroEigenVecs, tol=1.e-6)
        m1,m2 = tstools.minima_from_ts(pot.getEnergyGradient, ret.coords, ret.eigenvec, displace=1e-2)
        print "Energies: ", m1[1],ret.energy,m2[1]
        return [ret.coords,ret.energy],m1,m2
    
    def create_local_connect(self):
        mindist = self.create_mindist_object()
        pot = self.create_potential()
        
        local_connect = LocalConnect(pot, mindist)
        return local_connect
    
    def smooth_path(self, path, **kwargs):
        mindist = self.create_mindist_object()
        return smoothPath(path, mindist, **kwargs)
    
    def create_double_ended_connect(self, min1, min2, database):
        mindist = self.create_mindist_object()
        pot = self.create_potential()
        return DoubleEndedConnect(min1, min2, pot, mindist, database)


        

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
