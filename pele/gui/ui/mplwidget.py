from PyQt4 import QtGui, QtCore
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from PyQt4.QtGui import QWidget, QVBoxLayout

from matplotlib.figure import Figure
import numpy as np

from pele.utils.disconnectivity_graph import DisconnectivityGraph
from pele.storage import Database


class MPLWidget(FigureCanvas):
    """ defines a matplotlib widget """
    def __init__(self, parent=None): #, width=5, height=4, dpi=100):
        self.create_figure()
#        self.compute_initial_figure()
        
        FigureCanvas.__init__(self, self.fig)
        
        if parent is not None:
            self.setParent(parent)
        #self.reparent(parent, QtCore.QPoint(0, 0))

#        FigureCanvas.setSizePolicy(self,
#                                   QtGui.QSizePolicy.Expanding,
#                                   QtGui.QSizePolicy.Expanding)
#        FigureCanvas.updateGeometry(self)

#    def compute_initial_figure(self):
#        t = np.arange(0.0, 3.0, 0.01)
#        s = np.sin(2*np.pi*t)
#        self.axes.plot(t, s)


    def create_figure(self):
        self.fig = Figure(facecolor="white") #figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        self.axes.hold(True)
        



#    def sizeHint(self):
#        w, h = self.get_width_height()
#        return QtCore.QSize(0,0)
#        return QtCore.QSize(w, h)

#    def minimumSizeHint(self):
#        return QtCore.QSize(10, 10)

class MPLWidgetWithToolbar(QWidget):
    """ defines a matplotlib widget """
    def __init__(self, *args, **kwargs):
        super(MPLWidgetWithToolbar, self).__init__(*args, **kwargs)
        
        self.canvas = MPLWidget(parent=self)
        
        self.mpl_toolbar = NavigationToolbar(self.canvas, self)
        
        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.mpl_toolbar)
        self.setLayout(vbox)
        

        
        
        
