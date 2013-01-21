import numpy as np
from PyQt4.QtGui import QDialog, QApplication
import dgraph_browser
import sys

from pygmin.utils.disconnectivity_graph import DisconnectivityGraph
from pygmin.landscape import Graph
from pygmin.storage import Database
from pygmin.utils.events import Signal

class DGraphDialog(QDialog):
    def __init__(self, database):
        super(DGraphDialog, self).__init__()
        
        self.database = database
        
        self.ui = dgraph_browser.Ui_Form()
        self.ui.setupUi(self)
        self.plw = self.ui.widget
            
        self.params = {}
        self.set_defaults()
        self.rebuild_disconnectivity_graph()
        
        self.minimum_selected = Signal()
        # self.minimum_selected(minim)

    def set_defaults(self):
        self.ui.chkbx_center_gmin.setChecked(True)
        self.ui.chkbx_show_minima.setChecked(True)
        self.ui.chkbx_order_by_energy.setChecked(False)
        self.ui.chkbx_order_by_basin_size.setChecked(True)
        self.ui.chkbx_include_gmin.setChecked(True)


    def _get_input_parameters(self):
        self.params = {}
        params = self.params
        
        Emax = self.ui.lineEdit_Emax.text()
        if len(Emax) > 0:
            self.params["Emax"] = float(Emax)

        subgraph_size = self.ui.lineEdit_subgraph_size.text()
        if len(subgraph_size) > 0:
            self.params["subgraph_size"] = int(subgraph_size)

        nlevels = self.ui.lineEdit_nlevels.text()
        if len(nlevels) > 0:
            self.params["nlevels"] = int(nlevels)


        params["center_gmin"] = self.ui.chkbx_center_gmin.isChecked()
        params["show_minima"] = self.ui.chkbx_show_minima.isChecked()
        params["order_by_energy"] = self.ui.chkbx_order_by_energy.isChecked()
        params["order_by_basin_size"] = self.ui.chkbx_order_by_basin_size.isChecked()
        params["include_gmin"] = self.ui.chkbx_include_gmin.isChecked()


    def rebuild_disconnectivity_graph(self):
        
        self._get_input_parameters()
        self._build_disconnectivity_graph(**self.params)


    def _build_disconnectivity_graph(self, show_minima=True, **params):
        #this should be somewhere else
        db = self.database
        ax = self.plw.axes
        
        ax.clear()
        ax.hold(True)
        
        #make the plot look prettier
        ax.tick_params(axis='y', direction='out')
        ax.yaxis.tick_left()
        ax.spines['left'].set_color('black')
        ax.spines['left'].set_linewidth(0.5)
        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_color('none')
        ax.spines['right'].set_color('none')

        
        
        
        graphwrapper = Graph(db)
        dg = DisconnectivityGraph(graphwrapper.graph, **params)
        dg.calculate()
        
        
        #draw minima as points
        if show_minima:
            xpos, minima = dg.get_minima_layout()
            energies = [m.energy for m in minima]
            points = ax.scatter(xpos, energies, picker=5)
        
        
        def on_pick(event):
            if event.artist != points:
#                print "you clicked on something other than a node"
                return True
            thispoint = event.artist
            ind = event.ind[0]
            min1 = minima[ind]
            print "you clicked on minimum with id", min1._id, "and energy", min1.energy
            self.minimum_selected(min1)
        self.plw.fig.canvas.mpl_connect('pick_event', on_pick)

        
        
        
        #draw line segments connecting minima
        line_segments = dg.line_segments
        for x, y in line_segments:
            ax.plot(x, y, 'k')
        
        #adjust the limits of the graph.  this might be slow
        ymax = max(max(y) for x, y in line_segments)
        ymin = min(min(y) for x, y in line_segments)
        xmax = max(max(x) for x, y in line_segments)
        xmin = min(min(x) for x, y in line_segments)
        d = .2
        xmax += d
        xmin -= d
        ymax += d
        ymin -= d*4.
        self.plw.axes.set_xlim(xmin, xmax)
        self.plw.axes.set_ylim(ymin, ymax)
        
        
        ax.set_xticks([])

        self.plw.draw()



if __name__ == "__main__":
    
    db = Database("test.db")
    
    app = QApplication(sys.argv)        
    md = DGraphDialog(db)
    md.show()
    
    sys.exit(app.exec_()) 
        