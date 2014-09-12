from PyQt4 import QtGui
from ui.ui_takestep_explorer import Ui_MainWindow as UI
import numpy as np

class QMinimumInList(QtGui.QListWidgetItem):
    
    def __init__(self, minimum):
        text="%.4f (%d)"%(minimum.energy, minimum._id)
        QtGui.QListWidgetItem.__init__(self, text)
        self.minimum = minimum
        
    def __lt__(self, item2):
        #sort the energies in the list lowest to highest
        return self.minimum.energy > item2.minimum.energy

class TakestepExplorer(QtGui.QMainWindow):
    def __init__(self, parent=None, system=None, app=None, database=None):
        QtGui.QMainWindow.__init__(self, parent=parent)
    
        self.ui = UI()
        self.ui.setupUi(self)
        
        self.system = system
        self.app = app
        self.database = database
        self.coords = None        
        self.quenched = None
        
        self.ui.show3d.setSystem(system)
        
        self.database.on_minimum_added.connect(self.NewMinimum)
        self.read_minima()
        self.update_view()

    def read_minima(self):
        for minimum in self.system.database.minima():
            self.NewMinimum(minimum, sort_items=False)
        self.ui.listMinima.sortItems(1)

    def NewMinimum(self, minimum, sort_items=True):
        """ add a new minimum to the system """
        item = QMinimumInList(minimum)
        self.ui.listMinima.addItem(item)    
        if sort_items:
            self.ui.listMinima.sortItems(1)
        
    def on_actionDisplace_triggered(self, checked=None):
        if  checked is None or self.quenched is None:
            return
        if self.coords is None:
            self.coords = self.quenched.copy()
            
        
        takestep = self.system.get_takestep()
        takestep.takeStep(self.coords)
        self.update_view()
    
    def quench_event(self, coords=None, **kwargs):
        self.quench_path.append(coords.copy())
        
    def on_actionQuench_triggered(self, checked=None):
        if  checked is None or self.coords is None:
            return
        
        get_path = self.ui.actionShow_path.isChecked()
        
        if get_path:
            events = [self.quench_event]
        else:
            events = []
        self.quench_path = []
        quencher = self.system.get_minimizer(events=events)
        ret = quencher(self.coords)
        coords = ret[0]
        E = ret[1]
        print "energy", E
        self.quenched = coords
        self.coords = None
        
        self.database.addMinimum(E, coords)
        self.update_view(with_path=get_path)
        
    def on_listMinima_currentItemChanged(self, new, old):
        self.coords = None
        self.quenched = new.minimum.coords
        self.update_view()
        
    def update_view(self, with_path=False):
        pot = self.system.get_potential()        
        label = ""
        if self.quenched is not None:
            e, grad = pot.getEnergyGradient(self.quenched)
            label = "quenched: energy = %f, rms = %s\n"%(e, np.linalg.norm(grad)/np.sqrt(grad.size)) 
        if self.coords is not None:
            e, grad = pot.getEnergyGradient(self.coords)
            label += "instant: energy = %f, rms = %s"%(e, np.linalg.norm(grad)/np.sqrt(grad.size)) 
        
        self.ui.label.setText(label)
        if with_path:
            coordspath = np.array(self.quench_path)
            self.ui.show3d.setCoordsPath(coordspath, frame=-1)
        else:
            self.ui.show3d.setCoords(self.coords, index=2)
            self.ui.show3d.setCoords(self.quenched, index=1)
        

if __name__ == "__main__":
    import sys
    import pylab as pl
    from OpenGL.GLUT import glutInit
    glutInit()
    app = QtGui.QApplication(sys.argv)
    from pele.systems import LJCluster
    pl.ion()
    natoms = 13
    system = LJCluster(natoms)
    system.params.double_ended_connect.local_connect_params.NEBparams.iter_density = 5.
    x1, e1 = system.get_random_minimized_configuration()[:2]
    x2, e2 = system.get_random_minimized_configuration()[:2]
    db = system.create_database(db=":memory:")
    system.database = db
    min1 = db.addMinimum(e1, x1)
    min2 = db.addMinimum(e2, x2)
    
    wnd = TakestepExplorer(app=app, system=system, database = db)
    wnd.show()
    from PyQt4.QtCore import QTimer
    sys.exit(app.exec_()) 