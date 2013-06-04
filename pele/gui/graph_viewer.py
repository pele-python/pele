import networkx as nx
import numpy as np
from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import QDialog, QWidget

from pele.gui.ui.mplwidget import MPLWidget
from pele.landscape import TSGraph
from pele.utils.events import Signal

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s


class GraphViewWidget(MPLWidget):
    def __init__(self, database=None, parent=None, app=None, minima=None):
        MPLWidget.__init__(self, parent=parent)
#        self.widget = GraphDisplayWidget(parent=parent)
        self.database = database
        self.minima = minima
        self.app = app
        
        self.on_minima_picked = Signal()
        
    def make_graph(self, database=None, minima=None):
        if database is None:
            database = self.database
        if minima is None:
            minima = self.minima

        print "making graph", database, minima
        #get the graph object, eliminate nodes without edges
        graphwrapper = TSGraph(database, minima)
        graph = graphwrapper.graph
        print graph.number_of_nodes()
        degree = graph.degree()
        nodes = [n for n, nedges in degree.items() if nedges > 0]
        self.graph = graph.subgraph(nodes)
        print self.graph.number_of_nodes(), self.graph.number_of_edges()
        
    
    def show_graph(self):
        import pylab as pl
        if not hasattr(self, "graph"):
            self.make_graph()
        
        print "showing graph"
        ax = self.axes
        ax.clear()
        graph = self.graph
        
        #get the layout of the nodes from networkx
        layout = nx.spring_layout(graph)
        layoutlist = layout.items()
        xypos = np.array([xy for n, xy in layoutlist])
        #color the nodes by energy
        e = np.array([m.energy for m, xy in layoutlist])
        #plot the nodes
        points = ax.scatter(xypos[:,0], xypos[:,1], picker=5, 
                            s=8**2, c=e, cmap=pl.cm.autumn)
        if not hasattr(self, "colorbar"):
            self.colorbar = self.fig.colorbar(points)
        
        # label the nodes
        ids = [n._id for n, xy in layoutlist]
        for i in range(len(ids)):
            ax.annotate( ids[i], xypos[i] )
        
        #plot the edges as lines
        for u, v in graph.edges():
            line = np.array([layout[u], layout[v]])
            ax.plot(line[:,0], line[:,1], '-k')
        
        #scale the axes so the points are not cutoff
        xmin, ymin = np.min(xypos, 0)
        xmax, ymax = np.max(xypos, 0)
        dx = (xmax - xmin)*.1
        dy = (ymax - ymin)*.1
        ax.set_xlim([xmin-dx, xmax+dx])
        ax.set_ylim([ymin-dy, ymax+dy])
        
        
        def on_pick(event):
            if event.artist != points:
#                print "you clicked on something other than a node"
                return True
            ind = event.ind[0]
            min1 = layoutlist[ind][0]
            print "you clicked on minimum with id", min1._id, "and energy", min1.energy
            self.on_minima_picked(min1)

        
        self.fig.canvas.mpl_connect('pick_event', on_pick)

        
        self.draw()
        self.app.processEvents()
        

class GraphViewDialog(QDialog):
    def __init__(self, database, parent=None, app=None):
        QDialog.__init__(self, parent=parent)

        self.widget = GraphViewWidget(database=database, parent=self, app=app)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        self.widget.setSizePolicy(sizePolicy)
        self.widget.setObjectName(_fromUtf8("widget"))
#        self.widget.show_graph()

        self.retranslateUi(self)
#        QtCore.QMetaObject.connectSlotsByName(self)
        
        self.app = app

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(QtGui.QApplication.translate("Dialog", "Dialog", None, QtGui.QApplication.UnicodeUTF8))

    def start(self):
        self.widget.show_graph()
        self.widget.show_graph()

def start():
    wnd.start()


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
    
    #get some minima
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
    
        
    
    
    wnd = GraphViewDialog(db, app=app)
#    decrunner = DECRunner(system, db, min1, min2, outstream=wnd.textEdit_writer)
    glutInit()
    wnd.show()
    from PyQt4.QtCore import QTimer
    QTimer.singleShot(10, start)
    sys.exit(app.exec_()) 
