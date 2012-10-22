from PyQt4 import QtGui
import NewLJ
import numpy as np
from pygmin.NEB import NEB,dimer,tstools,InterpolatedPath
import pygmin.potentials.lj as lj
from pygmin import gui
                
class LJSystem(gui.GUISystem):
    def __init__(self):
        dlg = NewLJDialog()
        dlg.exec_()
        self.natoms = dlg.natoms()
        if dlg.result() == QtGui.QDialog.Rejected:
            raise BaseException("Aborted parameter dialog")
        
    def create_basinhopping(self):
        import pygmin.basinhopping as bh
        from pygmin.takestep import displace
        coords = np.random.random(3 * self.natoms)
        potential = lj.LJ()
        step = displace.RandomDisplacement(stepsize=0.5)
        opt = bh.BasinHopping(coords,potential,
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
    
    def Align(self, coords1, coords2):
        from pygmin.mindist.minpermdist_stochastic import minPermDistStochastic as minpermdist
        dist, X1, X2 = minpermdist( coords1, coords2, niter = 100 )
        return X1, X2
    
    def createNEB(self, coords1, coords2):
        return NEB.NEB(InterpolatedPath(coords1, coords2, 30), lj.LJ(), k = 100.)
    
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
