from PyQt4 import QtGui, QtCore
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

class MPLWidget(FigureCanvas):
    ''' defines a matplotlib widget '''
    def __init__(self, parent=None): #, width=5, height=4, dpi=100):
        self.create_figure() 
        FigureCanvas.__init__(self, self.fig)
        
        #self.reparent(parent, QtCore.QPoint(0, 0))

        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def create_figure(self):
        self.fig = Figure(facecolor="white") #figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        self.axes.hold(False)


    def sizeHint(self):
        w, h = self.get_width_height()
        return QtCore.QSize(w, h)

    def minimumSizeHint(self):
        return QtCore.QSize(10, 10)