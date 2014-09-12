import numpy as np

from PyQt4 import QtGui, QtCore

from pele.gui.ui.cv_viewer_ui import Ui_Form
from pele.thermodynamics import GetThermodynamicInfoParallel, minima_to_cv
from pele.utils.events import Signal

class GetThermodynamicInfoParallelQT(GetThermodynamicInfoParallel):
    def __init__(self, *args, **kwargs):
        super(GetThermodynamicInfoParallelQT, self).__init__(*args, **kwargs)
        self.on_finish = Signal()

    def poll(self):
        if self.njobs == 0:
            self.refresh_timer.stop()
            self.finish()
            return
        if not self.done_queue.empty():
            self.njobs -= 1
            ret = self.done_queue.get()
            self._process_return_value(ret)
    
    def finish(self):
        super(GetThermodynamicInfoParallelQT, self).finish()
        self.on_finish()
     
    def start(self):
        # populate the queue
        self._populate_queue()

        # start the workers
        for worker in self.workers:
            worker.start()
        
        self.refresh_timer = QtCore.QTimer()
        self.refresh_timer.timeout.connect(self.poll)
        self.refresh_timer.start(50.) # time in msec

class HeatCapacityWidget(QtGui.QWidget):
    def __init__(self, system, database, parent=None):
        super(HeatCapacityWidget, self).__init__(parent=parent)
        self.ui = Ui_Form()
        self.ui.setupUi(self)
        
        self.system = system
        self.database = database
        
        self.canvas = self.ui.mplwidget.canvas
        self.axes = self.canvas.axes
    
    def rebuild_cv_plot(self):
        self._compute_thermodynamic_info(on_finish=self.make_cv_plot)
    
    def make_cv_plot(self):
        self._compute_cv()
        self._plot_cv()
    
    def _get_ndof(self):
        return self.system.get_ndof()
        
    def _get_nmin_max(self):
        txt = self.ui.lineEdit_nmin_max.text()
        if len(txt) > 0:
            nmin_max = int(txt)
        else:
            nmin_max = None
        return nmin_max

    def _compute_thermodynamic_info(self, nproc=2, on_finish=None, verbose=False):
        nmin = self._get_nmin_max()
        if nmin is not None and nmin < self.database.number_of_minima():
            self.minima = self.database.minima()[:nmin]
        else:
            self.minima = self.database.minima()

        self.worker = GetThermodynamicInfoParallelQT(self.system, self.database, 
                            npar=nproc, verbose=verbose, only_minima=True)
        if on_finish is not None:
            self.worker.on_finish.connect(on_finish)
        self.worker.start()
        
        njobs = self.worker.njobs
        self.ui.label_status.setText("computing thermodynamic information for %d minima" % njobs)
    
    
    def _compute_cv(self):
        Tlist = self._get_T_range()
        lZ, U, U2, Cv = minima_to_cv(self.minima, Tlist, self._get_ndof())
        self.Tlist = Tlist
        self.Cv = Cv
    
    def _get_Tmin(self):
        res = 0.01        
        txt = self.ui.lineEdit_Tmin.text()
        if len(txt) > 0:
            res = float(txt)
        return res
    
    def _get_Tmax(self):
        res = 1.        
        txt = self.ui.lineEdit_Tmax.text()
        if len(txt) > 0:
            res = float(txt)
        return res
    def _get_nT(self):
        res = 100        
        txt = self.ui.lineEdit_nT.text()
        if len(txt) > 0:
            res = float(txt)
        return res

    
    def _get_T_range(self):
        Tmin = self._get_Tmin()
        Tmax = self._get_Tmax()
        nT = self._get_nT()
        
        dT = (Tmax - Tmin) / nT
        return np.arange(Tmin, Tmax, dT)
    
    def _plot_cv(self):
        self.ui.label_status.setText("showing heat capacity calculated with %d minima" % len(self.minima))
        axes = self.axes
        axes.clear()
        axes.plot(self.Tlist, self.Cv)
        self.canvas.draw()
    
    def on_btn_recalculate_clicked(self, clicked=None):
        if clicked is None: return 
        self.rebuild_cv_plot()

class HeatCapacityViewer(QtGui.QMainWindow):
    def __init__(self, system, database, parent=None, app=None):
        super(HeatCapacityViewer, self).__init__( parent=parent)
        self.cv_widget = HeatCapacityWidget(system, database, parent=self)
        self.setCentralWidget(self.cv_widget)
        self.setWindowTitle("Harmonic Superposition Heat Capacity")
    
    def rebuild_cv_plot(self):
        self.cv_widget.rebuild_cv_plot()
        



def test():
    import sys
    from pele.systems import LJCluster
    app = QtGui.QApplication(sys.argv)
    system = LJCluster(13)
    
    db = system.create_database()
    bh = system.get_basinhopping(db, outstream=None)
    bh.run(200)

    obj = HeatCapacityViewer(system, db)
    obj.show()

    def test_start():
        obj.rebuild_cv_plot()
    
    QtCore.QTimer.singleShot(10, test_start)
    sys.exit(app.exec_()) 




if __name__ == "__main__":
    test()
        