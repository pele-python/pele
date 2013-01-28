from PyQt4 import QtGui, QtCore, Qt
from PyQt4.QtGui import QDockWidget
from pygmin.utils.events import Signal

from nebdlg import NEBWidget
from show3d import Show3D
from pygmin.storage import Database
from pygmin.gui.ui.mplwidget import MPLWidget
from  pygmin.gui.ui.ui_neb_explorer import Ui_MainWindow as UI
from dlg_params import DlgParams

class NEBRunner(object):
    def __init__(self, app, system, freq = 30):
        self.system = system
        self.on_update_gui = Signal()
        self.frq = freq
        self.app = app
    
    def run(self, coords1, coords2):
        self.count=0
        neb = self.create_neb(coords1, coords2)
        self.neb = neb
        neb.update_event.connect(self._neb_update)
        
        self.k = []
        self.nimages = []
        self.energies=[]
        self.stepnum = []
        self.distances = []
        self.rms = []
        neb.run()
        
    def _neb_update(self, energies=None, distances=None, stepnum=None, path=None, rms=None, **kwargs):
        self.app.processEvents()
        self.count += 1
        if self.count % self.frq == 1:
            self.stepnum.append(stepnum)
            self.k.append(self.neb.neb.k)
            self.rms.append(rms)
            self.nimages.append(len(self.neb.neb.coords))
            self.energies.append(energies.copy())
            self.distances.append(distances.copy())
            self.path = path.copy()
            self.on_update_gui(self)
                       
    def create_neb(self, coords1, coords2):
        """setup the NEB object"""
        system = self.system
        
        throwaway_db = Database()
        min1 = throwaway_db.addMinimum(0., coords1)
        min2 = throwaway_db.addMinimum(1., coords2)
        #use the functions in DoubleEndedConnect to set up the NEB in the proper way
        double_ended = system.get_double_ended_connect(min1, min2, 
                                                            throwaway_db, 
                                                            fresh_connect=True)
        local_connect = double_ended._getLocalConnectObject()
    
        
        
        return local_connect.create_neb(system.get_potential(),
                                          coords1, coords2,
                                          **local_connect.NEBparams)        

class NEBEnergyWidget(MPLWidget):
    def __init__(self, parent=None, nplots = 3):
        #QtGui.QWidget
        MPLWidget.__init__(self, parent=parent)
        #self.canvas = MPLWidget(self)
        self.nplots = nplots
        
    def update_gui(self, nebrunner):
        nplots = min(self.nplots, len(nebrunner.energies))
        self.axes.clear()
        self.axes.set_xlabel("distance")
        self.axes.set_ylabel("energy")
        
        for distances, energies, stepnum in zip(nebrunner.distances[-nplots:],
                                       nebrunner.energies[-nplots:], nebrunner.stepnum[-nplots:]):
            acc_dist = [0.]
            acc=0.
            for d in distances:
                acc+=d
                acc_dist.append(acc)
            self.axes.plot(acc_dist, energies, "o-", label="step %d"%stepnum)
        self.axes.legend(loc='best')
        self.draw()

class NEBDistanceWidget(MPLWidget):
    def __init__(self, parent=None):
        #QtGui.QWidget
        MPLWidget.__init__(self, parent=parent)
        #self.canvas = MPLWidget(self)
        
    def update_gui(self, nebrunner):
        self.axes.clear()
        self.axes.set_xlabel("relative distance")
        self.axes.set_ylabel("image")
        self.axes.set_title("distances")
        self.axes.plot(nebrunner.distances[-1])
        self.draw()
        
class NEBTimeseries(MPLWidget):
    def __init__(self, parent=None, attrname="k", yscale='linear'):
        MPLWidget.__init__(self, parent=parent)
        self.neb_attribute=attrname
        self.yscale = yscale
        
    def update_gui(self, nebrunner):
        value = getattr(nebrunner, self.neb_attribute)
        stepnum = nebrunner.stepnum
        
        self.axes.clear()
        self.axes.set_xlabel("step")
        self.axes.set_ylabel(self.neb_attribute)
        self.axes.set_yscale(self.yscale)
        self.axes.plot(stepnum, value)
        self.draw()
        
class NEBExplorer(QtGui.QMainWindow):
    def __init__(self, parent=None, system=None, app=None):
        QtGui.QMainWindow.__init__(self)
    
        self.ui = UI()
        self.ui.setupUi(self)
        
        self.system = system
        self.mdi = QtGui.QMdiArea(self)
        self.setCentralWidget(self.mdi)
        
        self.nebrunner = NEBRunner(app, system)
#        
#        from dlg_params import EditParamsWidget
#        w = QtGui.QDockWidget("NEB parameters", self)
#        w.setWidget(EditParamsWidget(self, 
#                         self.system.params.double_ended_connect.local_connect_params.NEBparams))
#        self.addDockWidget(QtCore.Qt.RightDockWidgetArea, w)
#        
#        self.editparams = w
#        
        
        self.view_energies = self.new_view("Energies", NEBEnergyWidget(), QtCore.Qt.TopDockWidgetArea)
        self.view_distances = self.new_view("Distances", NEBDistanceWidget(), QtCore.Qt.TopDockWidgetArea)
        self.view_k = self.new_view("k", NEBTimeseries(attrname="k"), QtCore.Qt.BottomDockWidgetArea)
        self.view_nimages = self.new_view("nimages", NEBTimeseries(attrname="nimages"), QtCore.Qt.BottomDockWidgetArea)
        self.view_rms = self.new_view("rms", NEBTimeseries(attrname="rms", yscale='log'), QtCore.Qt.BottomDockWidgetArea)
        
        
        self.centralWidget().hide()
                    
    def new_view(self, title, widget, pos=QtCore.Qt.RightDockWidgetArea):
        child = QtGui.QDockWidget(title, self)
        child.setWidget(widget)
        self.addDockWidget(pos, child)
        self.nebrunner.on_update_gui.connect(widget.update_gui)
        return child
    
    def new_neb(self, coords1, coords2):
        self.coords1 = coords1.copy()
        self.coords2 = coords2.copy()
        self.nebrunner.run(coords1, coords2)
        
    def toggle_view(self, view, show):
        if show:
            view.show()
        else:
            view.hide()

    def on_actionRun_triggered(self):
        self.nebrunner.run(self.coords1, self.coords2)
        
    def on_actionParams_triggered(self):
        if not hasattr(self, "paramsdlg"):
            self.paramsdlg = DlgParams(self.system.params.double_ended_connect.local_connect_params.NEBparams)
        self.paramsdlg.show()
        
    def on_actionRms_toggled(self, checked):
        self.toggle_view(self.view_rms, checked)
    def on_actionE_toggled(self, checked):
        self.toggle_view(self.view_energies, checked)
    def on_actionS_toggled(self, checked):
        self.toggle_view(self.view_distances, checked)
    def on_actionK_toggled(self, checked):
        self.toggle_view(self.view_k, checked)
    def on_actionNimages_toggled(self, checked):
        self.toggle_view(self.view_nimages, checked)
        
def start():
    wnd.new_neb(x1, x2)
    
if __name__ == "__main__":
    import sys
    import pylab as pl
    
    app = QtGui.QApplication(sys.argv)
    from pygmin.systems import LJCluster
    pl.ion()
    natoms = 13
    system = LJCluster(natoms)
    system.params.double_ended_connect.local_connect_params.NEBparams.iter_density = 5.
    x1, e1 = system.get_random_minimized_configuration()[:2]
    x2, e2 = system.get_random_minimized_configuration()[:2]
    db = Database()
    min1 = db.addMinimum(e1, x1)
    min2 = db.addMinimum(e2, x2)
    
    wnd = NEBExplorer(app=app, system=system)
    wnd.show()
    from PyQt4.QtCore import QTimer
    QTimer.singleShot(10, start)
    sys.exit(app.exec_()) 