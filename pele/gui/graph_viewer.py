import networkx as nx
import numpy as np
from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import QDialog, QWidget

#from pele.gui.ui.mplwidget import MPLWidgetWithToolbar
from pele.gui.ui.graph_view_ui import Ui_Form
from pele.landscape import TSGraph
from pele.utils.events import Signal

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s


class GraphViewWidget(QWidget):
    def __init__(self, database=None, parent=None, app=None, minima=None):
        QWidget.__init__(self, parent=parent)
        self.ui = Ui_Form()
        self.ui.setupUi(self)
#        self.widget = GraphDisplayWidget(parent=parent)
        self.database = database
        self.minima = minima
        self.app = app
        
        self.canvas = self.ui.canvas.canvas
        
        self.axes = self.canvas.axes
        self.fig = self.canvas.fig
        
        self.on_minima_picked = Signal()
        
        self.from_minima = set()
        self.positions = dict()
    
    
    
    def on_btn_show_all_clicked(self, clicked=None):
        if clicked is None: return
        self.show_all()
        
    def show_all(self):
        self.ui.label_status.setText("showing full graph")
        self.from_minima.clear()
        self.positions.clear()
        self.make_graph()
        self.show_graph()
    
    def make_graph_from(self, minima, cutoff=1):
        cutoff = cutoff + 1
        self.from_minima.update(minima)
        minima = self.from_minima
        graph = self.full_graph
        nodes = set()
        # make a graph from the minima in self.minima and nearest neighbors
        for m in minima:
            nodesdir = nx.single_source_shortest_path(graph, m, cutoff=cutoff)
            nodes.update(nodesdir)
        
        
        self.graph = graph.subgraph(nodes)

        # remove nodes not in the graph from the dictionary positions
        difference = set(self.positions.viewkeys())
        difference.difference_update(self.graph.nodes())
        for m in difference:
            self.positions.pop(m)
    
    def make_graph(self, database=None, minima=None):
        if database is None:
            database = self.database
        if minima is None:
            minima = self.minima

        print "making graph", database, minima
        #get the graph object, eliminate nodes without edges
        graphwrapper = TSGraph(database, minima)
        graph = graphwrapper.graph
        self.full_graph = graph
        print graph.number_of_nodes()
        degree = graph.degree()
        nodes = [n for n, nedges in degree.items() if nedges > 0]
        self.graph = graph.subgraph(nodes)
        print self.graph.number_of_nodes(), self.graph.number_of_edges()
    
    def on_minima_picked_local(self, min1):
        if self.ui.checkBox_zoom.isChecked():
            self.make_graph_from([min1])
            text = "showing graph near minima "
            for m in self.from_minima: 
                text += " " + str(m._id)
            self.ui.label_status.setText(text)
            self.show_graph()
    
    def show_graph(self, fixed=False):
        import pylab as pl
        if not hasattr(self, "graph"):
            self.make_graph()
        
        print "showing graph"
        ax = self.axes
        ax.clear()
        graph = self.graph
        
        #get the layout of the nodes from networkx
        oldlayout = self.positions
        layout = nx.spring_layout(graph, pos=oldlayout)#, fixed=fixed)
        self.positions.update(layout)
        layout = self.positions

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
            self.on_minima_picked_local(min1)

        
        self.fig.canvas.mpl_connect('pick_event', on_pick)

        
        self.canvas.draw()
        self.app.processEvents()
        

class GraphViewDialog(QtGui.QMainWindow):
    def __init__(self, database, parent=None, app=None):
        QtGui.QMainWindow.__init__(self, parent=parent)

        self.widget = GraphViewWidget(database=database, parent=self, app=app)
        self.setCentralWidget(self.widget)
        
        self.app = app


    def start(self):
        self.widget.show_all()
#        self.widget.make_graph()
#        gmin = self.widget.database.minima()[0:1]
#        self.widget.make_graph_from(gmin)
#        self.widget.show_graph()

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
