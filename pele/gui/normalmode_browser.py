import pickle

import numpy as np
from PyQt4 import QtGui

from pele.gui.ui.ui_normalmode_explorer import Ui_MainWindow as UI
from pele.thermodynamics import normalmodes 
from pele.gui.dlg_params import DlgParams

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
    """
    the GUI for exploring normal modes
    """
    def __init__(self, parent=None, system=None, app=None):
        QtGui.QMainWindow.__init__(self, parent=parent)
        
        self.ui = UI()
        self.ui.setupUi(self)
        
        self.ui.view3D.setSystem(system)
        self.system = system
        self._params = dict()
        self._params["amplitude"]=1.0
        self._params["remove_known_zeroev"]=True
        export = self._params["export"] = dict()
        export["nframes"]=100
        
        self._params["nframes"] = 30
        
        self.app = app
        self.current_selection = None
        
        self.ui.actionShow_energies.setChecked(False)
        self.ui.mplwidget.hide()
        
        self.ui.actionRun.setVisible(False)
         
    def set_coords(self, coords, normalmodes=None):
        """
        set the coordinates for which the normal modes will be computed
        """
        self.coords = coords
        self.normalmodes = normalmodes        
        
        if normalmodes is None:
            self._calculate_normalmodes()
        
        self._fill_normalmodes()
        
        self.ui.view3D.setCoords(coords)
        
        if self.ui.actionShow_energies.isChecked():
            self.draw_energy_plot()


    def _use_hessian_eigs(self):
        return self.ui.actionHessian_eigs.isChecked()
    
    def on_actionHessian_eigs_toggled(self, checked):
        if checked is None: return
        self.set_coords(self.coords)
        
    
    def _calculate_normalmodes(self):
        """
        compute the normal modes
        """
        if self._use_hessian_eigs():
            pot = self.system.get_potential()
            hess = pot.getHessian(self.coords)
            freq, mode = normalmodes(hess, metric=None)
        else:
            freq, mode = self.system.get_normalmodes(self.coords)
        mode=np.real(mode.transpose())
        
        self.normalmodes = []
        #self.normalmodes.append((fre[0], m.flatten()))
        for f, m in zip(freq, mode):
            self.normalmodes.append((f, m)) #np.dot(metric, m)))
             
    def _fill_normalmodes(self):
        """
        populate the list of normal modes
        """
        row = self.ui.listNormalmodes.currentRow()
        self.ui.listNormalmodes.clear()
        for n in self.normalmodes:
            self.ui.listNormalmodes.addItem(NormalmodeItem(n))
        self.current_selection = self.ui.listNormalmodes.item(row)
        self.ui.listNormalmodes.setCurrentRow(row)

    def on_listNormalmodes_currentItemChanged(self, newsel):
        """
        change which normal mode we're looking at
        """
        if newsel is None:
            self.currentmode = None
            return
        orthogopt = self.system.get_orthogonalize_to_zero_eigenvectors()
        mode = newsel.get_mode().copy()
        if self._params["remove_known_zeroev"] and orthogopt is not None:
            mode = orthogopt(mode, self.coords)
         
        self.currentmode = mode
        self.current_selection = newsel
        
        # generate the configurations from the normal mode
        amp = self._params["amplitude"]
        vector = self.currentmode
        nframes = self._params["nframes"]
        dxlist = [amp * float(i) / nframes for i in xrange(-nframes/2,nframes/2)]
        coordspath = [self.coords + dx * vector for dx in dxlist] 
        coordspath = np.array(coordspath)
        
        self.dxlist = dxlist
        self.coordspath = coordspath
        
        self.ui.view3D.setCoordsPath(coordspath)#, labels=labels)
#        self.ui.view3D.ui.btn_animate.hide()

        if self.ui.actionShow_energies.isChecked():
            self.draw_energy_plot()
    
    def draw_energy_plot(self):
        """
        make a plot of the energies and the energies from the harmonic approximation
        """
        if self.current_selection is None:
            print "clearing energy axes"
            self.ui.mplwidget.axes.clear() 
            self.ui.mplwidget.draw()
            return
        dxlist = self.dxlist
        coordspath = self.coordspath
        
        # get the energies of the configurations
        pot = self.system.get_potential()
        energies = [pot.getEnergy(coords) for coords in coordspath]
        
        # get the energies of the harmonic approximation
        freq = self.current_selection.get_freq()
        expected_energies = np.array([ 0.5* freq * dx **2 for dx in dxlist])
        expected_energies += pot.getEnergy(self.coords)
        
        # make the plot
        ax = self.ui.mplwidget.axes
        ax.clear()
        ax.plot(dxlist, energies, label="energy")
        ax.plot(dxlist, expected_energies, label="harmonic approximation")
        ax.legend(loc='best')
        ax.set_xlabel("displacement")
        self.ui.mplwidget.draw()
    
    def on_actionRun_toggled(self, checked=None):
        if checked is None: return
        if checked:
            self.ui.view3D.start_animation()
        else:
            self.ui.view3D.stop_animation()
    
    def on_actionShow_energies_toggled(self, checked=None):
        if checked is None: return
        if checked:
            self.ui.mplwidget.show()
            self.draw_energy_plot()
        else:
            self.ui.mplwidget.hide()

    def on_actionSave_triggered(self, checked=None):
        """
        save the normal modes to disk
        """
        if checked is None:
            return
        dialog = QtGui.QFileDialog(self)
        dialog.setFileMode(QtGui.QFileDialog.AnyFile)
        dialog.selectFile("mode.pickle")
        dialog.setAcceptMode(QtGui.QFileDialog.AcceptSave)
        if not dialog.exec_():
            return
        filename = dialog.selectedFiles()[0]
        path = []
        nframes = self._params["export"]["nframes"]
        for i in xrange(nframes):
            t = np.sin(i/float(nframes)*2.*np.pi)
            path.append(self.coords + self._params["amplitude"]*t*self.currentmode)
        pickle.dump(path, open(filename, "w"))

    def on_actionParameters_triggered(self, checked=None):
        """
        open a dialog box to change the parameters
        """
        if checked is None:
            return
        if not hasattr(self, "_paramsdlg"):
            self._paramsdlg = DlgParams(self._params, parent=self)
        self._paramsdlg.show()            
    
        
if __name__ == "__main__":
    from OpenGL.GLUT import glutInit
    import sys
    glutInit()
    app = QtGui.QApplication(sys.argv)
    from pele.systems import LJCluster
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
