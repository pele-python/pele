from PyQt4 import QtGui

from pele.gui.ui.rate_gui import Ui_Form
from pele.rates import RateCalculation
from _cv_viewer import GetThermodynamicInfoParallelQT

class RateWidget(QtGui.QWidget):
    def __init__(self, system, database, temperature=1., parent=None):
        QtGui.QWidget.__init__(self, parent=parent)
        self.ui = Ui_Form()
        self.ui.setupUi(self)
        
        self.system = system
        self.database = database
        
        self.ui.lineEdit_T.setText(str(temperature))
        self.ui.lineEdit_A.setReadOnly(True)
        self.ui.lineEdit_B.setReadOnly(True)
        
        self.A = set()
        self.B = set()
    
    def update_selected_text(self):
        if len(self.A) > 0:
            m = iter(self.A).next()
            self.ui.lineEdit_A.setText("%s (%s)" % (m.energy, m._id))
        if len(self.B) > 0:
            m = iter(self.B).next()
            self.ui.lineEdit_B.setText("%s (%s)" % (m.energy, m._id))
        
    def update_A(self, minimum):
        self.A = {minimum}
        self.update_selected_text()
    def update_B(self, minimum):
        self.B = {minimum}
        self.update_selected_text()
    
#    def read_A(self):
#        line = self.ui.lineEdit_A.text()
#        A = set(map(int, line.split()))
#        return A
#    def read_B(self):
#        line = self.ui.lineEdit_B.text()
#        B = set(map(int, line.split()))
#        return B
    
    def _add_result(self, A, B, rAB, rBA):
        Aid = [m._id for m in A]
        Bid = [m._id for m in B]
        self.ui.textBrowser.append("rate %s -> %s = %s" % (Aid, Bid, rAB))
        self.ui.textBrowser.append("rate %s -> %s = %s" % (Bid, Aid, rBA))
        self.ui.textBrowser.append("")
    
    def _compute_rates(self):
        
        self.ui.label_status.setText("computing rates %s <-> %s" % (self.A, self.B))
        T = float(self.ui.lineEdit_T.text())
        calculator = RateCalculation(self.transition_states, self.A, self.B, 
                                     T=T, use_fvib=True)
        calculator.compute_rates()
        rAB = calculator.get_rate_AB()
        rBA = calculator.get_rate_BA()
        self._add_result(self.A, self.B, rAB, rBA)
        
        self.ui.label_status.setText("")
        
    
    def compute_rates(self):
        self._compute_thermodynamic_info(on_finish=self._compute_rates)
        
        

    def _compute_thermodynamic_info(self, nproc=2, on_finish=None, verbose=False):
        self.transition_states = list(self.database.transition_states())

        self.worker = GetThermodynamicInfoParallelQT(self.system, self.database, 
                            npar=nproc, verbose=verbose, only_minima=False)
        if on_finish is not None:
            self.worker.on_finish.connect(on_finish)
        self.worker.start()
        
        njobs = self.worker.njobs
        self.ui.label_status.setText("computing thermodynamic information for %d minima" % njobs)

    def on_btn_compute_clicked(self, clicked=None):
        if clicked is None: return
        self.compute_rates()

class RateViewer(QtGui.QMainWindow):
    def __init__(self, system, database, parent=None, app=None):
        super(RateViewer, self).__init__(parent=parent)
        self.rate_widget = RateWidget(system, database, parent=self)
        self.setCentralWidget(self.rate_widget)
        self.setWindowTitle("Rates")
        
        self.update_A = self.rate_widget.update_A
        self.update_B = self.rate_widget.update_B
        self.compute_rates = self.rate_widget.compute_rates
    
#    def rebuild_cv_plot(self):
#        self.cv_widget.rebuild_cv_plot()


def start():
    wnd.compute_rates()

if __name__ == "__main__":
    from OpenGL.GLUT import glutInit
    import sys
    import pylab as pl

    app = QtGui.QApplication(sys.argv)
    from pele.systems import LJCluster
    pl.ion()
    natoms = 13
    system = LJCluster(natoms)
    system.params.double_ended_connect.local_connect_params.NEBparams.iter_density = 5.
    dbname = "lj%dtest.db" % (natoms,)
    db = system.create_database(dbname)
    
    # get some minima
    if False:
        bh = system.get_basinhopping(database=db)
        bh.run(10)
        minima = db.minima()
    else:
        x1, e1 = system.get_random_minimized_configuration()[:2]
        x2, e2 = system.get_random_minimized_configuration()[:2]
        min1 = db.addMinimum(e1, x1)
        min2 = db.addMinimum(e2, x2)
        minima = [min1, min2]

    
    # connect some of the minima
    nmax = min(3, len(minima))
    m1 = minima[0]
    for m2 in minima[1:nmax]:
        connect = system.get_double_ended_connect(m1, m2, db)
        connect.connect()
    
        
    
    
    wnd = RateViewer(system, db)
#    decrunner = DECRunner(system, db, min1, min2, outstream=wnd.textEdit_writer)
    glutInit()
    wnd.show()
    from PyQt4.QtCore import QTimer
    wnd.update_A(db.minima()[0])
    wnd.update_B(db.minima()[1])
    
    QTimer.singleShot(10, start)
    sys.exit(app.exec_()) 

