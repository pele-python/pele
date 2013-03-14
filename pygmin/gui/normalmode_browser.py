from PyQt4 import QtGui, QtCore, Qt
from  ui.ui_normalmode_explorer import Ui_MainWindow as UI
from pygmin.thermodynamics import normalmodes 
import numpy as np
import pickle

class NormalmodeItem(QtGui.QListWidgetItem):    
    def __init__(self, normalmode):
        text="%.5e"%normalmode[0]
        self.normalmode = normalmode
        
        QtGui.QListWidgetItem.__init__(self, text)
        
    def get_mode(self):
        return self.normalmode[1]
    
    def get_freq(self):
        return self.normalmode[0]
    
    def __lt__(self, item2):
        #sort the energies in the list lowest to highest
        return self.normalmode[0] < item2.normalmode[0]
    
class NormalmodeBrowser(QtGui.QMainWindow):
    def __init__(self, parent=None, system=None, app=None):
        QtGui.QMainWindow.__init__(self, parent=parent)
        
        self.ui = UI()
        self.ui.setupUi(self)
        
        self.ui.view3D.setSystem(system)
        self.system = system
        
    def set_coords(self, coords, normalmodes=None):
        self.coords = coords
        self.normalmodes = normalmodes        
        
        if normalmodes is None:
            self._calculate_normalmodes()
        
        self._fill_normalmodes()
        
        self.ui.view3D.setCoords(coords)
        
    def _calculate_normalmodes(self):
        pot = self.system.get_potential()
        E, g, hess = pot.getEnergyGradientHessian(self.coords)
        metric = self.system.get_metric_tensor(self.coords)
        freq, mode = normalmodes(hess, metric = metric)
        mode=np.real(mode.transpose())
        
        self.normalmodes = []
        #self.normalmodes.append((fre[0], m.flatten()))
        for f, m in zip(freq, mode):
            self.normalmodes.append((f, m)) #np.dot(metric, m)))
             
    def _fill_normalmodes(self):
        for n in self.normalmodes:
            self.ui.listNormalmodes.addItem(NormalmodeItem(n))
        
    def on_listNormalmodes_currentItemChanged(self, newsel):
        self.currentmode = newsel.get_mode()
        
    def on_sliderFrame_valueChanged(self, val):
        displace = self.currentmode
        self.ui.view3D.setCoords(self.coords + displace*val/self.ui.sliderFrame.maximum())
        
    def on_actionRun_toggled(self, checked):
        if checked:
            self.animate=True
            self._animate_dir = 1
            QtCore.QTimer.singleShot(0., self._next_frame)
        else:
            self.animate = False
            
    def on_actionSave_triggered(self, checked=None):
        if checked is None:
            return
        dialog = QtGui.QFileDialog(self)
        dialog.setFileMode(QtGui.QFileDialog.AnyFile)
        dialog.selectFile("mode.pickle")
        dialog.setAcceptMode(QtGui.QFileDialog.AcceptSave);
        if(not dialog.exec_()):
            return
        filename = dialog.selectedFiles()[0]
        path = []
        for i in xrange(100):
            t = np.sin(i/100.*2.*np.pi)
            path.append(self.coords + t*self.currentmode)
        pickle.dump(path, open(filename, "w"))

    
    def _next_frame(self):
        cur = self.ui.sliderFrame.value()
        if cur == self.ui.sliderFrame.maximum():
            self._animate_dir = -1
        elif cur ==  self.ui.sliderFrame.minimum():
            self._animate_dir = 1
        cur +=self._animate_dir
        self.ui.sliderFrame.setValue(cur)
        self.on_sliderFrame_valueChanged(cur)
        
        if self.animate:
            QtCore.QTimer.singleShot(0.05, self._next_frame)
        
        
if __name__ == "__main__":
    from OpenGL.GLUT import glutInit
    import sys
    glutInit()
    app = QtGui.QApplication(sys.argv)
    from pygmin.systems import LJCluster
    natoms = 13
    system = LJCluster(natoms)
    system.params.double_ended_connect.local_connect_params.NEBparams.iter_density = 5.
    x1, e1 = system.get_random_minimized_configuration()[:2]
    db = system.create_database()
    match = system.get_compare_exact()
    min1 = db.addMinimum(e1, x1)
    
    com = match.measure.get_com(x1)
    match.transform.translate(x1, -com)
    wnd = NormalmodeBrowser(app=app, system=system)
    wnd.set_coords(x1)
    
    wnd.show()
    sys.exit(app.exec_())     
