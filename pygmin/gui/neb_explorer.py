from PyQt4 import QtGui, QtCore, Qt
from PyQt4.QtGui import QDockWidget
from pygmin.utils.events import Signal

from nebdlg import NEBWidget
from show3d import Show3D
from pygmin.storage import Database
from pygmin.gui.ui.mplwidget import MPLWidget

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
        neb.run()
        
    def _neb_update(self, energies=None, distances=None, stepnum=None, **kwargs):
        self.app.processEvents()
        self.count += 1
        if self.count % self.frq == 1:
            self.k = self.neb.neb.k
            self.nimages = len(self.neb.neb.coords)
            self.on_update_gui(energies = energies, distances=distances, 
                               stepnum=stepnum, nebrunner=self, **kwargs)
                       
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
    def __init__(self, parent=None):
        #QtGui.QWidget
        MPLWidget.__init__(self, parent=parent)
        #self.canvas = MPLWidget(self)
        
    def update_gui(self, energies=None, distances=None, stepnum=None, **kwargs):
        acc_dist = [0.]
        acc=0.
        for d in distances:
            acc+=d
            acc_dist.append(acc)
        self.axes.clear()
        self.axes.set_xlabel("distance")
        self.axes.set_ylabel("energy")
        self.axes.plot(acc_dist, energies)
        self.draw()

class NEBDistanceWidget(MPLWidget):
    def __init__(self, parent=None):
        #QtGui.QWidget
        MPLWidget.__init__(self, parent=parent)
        #self.canvas = MPLWidget(self)
        
    def update_gui(self, energies=None, distances=None, stepnum=None, **kwargs):
        self.axes.clear()
        self.axes.set_xlabel("distance")
        self.axes.set_ylabel("energy")
        self.axes.set_title("distances")
        self.axes.plot(distances)
        self.draw()
        
class NEBTimeseries(MPLWidget):
    def __init__(self, parent=None, attrname="k"):
        MPLWidget.__init__(self, parent=parent)
        self.neb_attribute=attrname
        self.timeseries = []
        self.stepnum = []
        
    def update_gui(self, energies=None, distances=None, stepnum=None, nebrunner=None, **kwargs):
        value = getattr(nebrunner, self.neb_attribute)
        self.timeseries.append(value)
        self.stepnum.append(stepnum)
        self.axes.clear()
        self.axes.set_xlabel("step")
        self.axes.set_ylabel(self.neb_attribute)
        self.axes.plot(self.stepnum, self.timeseries)
        self.draw()
        
class NEBExplorer(QtGui.QMainWindow):
    def __init__(self, parent=None, system=None, app=None):
        QtGui.QMainWindow.__init__(self)
    
        self.mdi = QtGui.QMdiArea(self)
        self.setCentralWidget(self.mdi)
        
        self.nebrunner = NEBRunner(app, system)
        
        self.new_view("Energies", NEBEnergyWidget(), QtCore.Qt.TopDockWidgetArea)
        self.new_view("Distances", NEBDistanceWidget(), QtCore.Qt.TopDockWidgetArea)
        self.new_view("k", NEBTimeseries(attrname="k"), QtCore.Qt.BottomDockWidgetArea)
        self.new_view("nimages", NEBTimeseries(attrname="nimages"), QtCore.Qt.BottomDockWidgetArea)
        #self.new_view("Distances", NEBEnergyWidget())
        #self.mdi.tileSubWindows()
        #self.new_doc("3d view", Show3D(self))
        self.centralWidget().hide()
                
#    def new_view(self, title, widget):
#        child = QtGui.QMdiSubWindow(self)
#        child.setWindowTitle(title)
#        self.mdi.addSubWindow(child)
#        child.setWidget(widget)
#        self.nebrunner.on_update_gui.connect(widget.update_gui)
#        return child
    
    def new_view(self, title, widget, pos=QtCore.Qt.RightDockWidgetArea):
        child = QtGui.QDockWidget(title, self)
        child.setWidget(widget)
        self.addDockWidget(pos, child)
        self.nebrunner.on_update_gui.connect(widget.update_gui)
        return child
    
    def new_neb(self, coords1, coords2):
        self.nebrunner.run(coords1, coords2)

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