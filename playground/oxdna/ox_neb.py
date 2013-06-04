import sys
from OpenGL.GLUT import glutInit
from oxgui import OXDNASystem
from PyQt4.QtGui import QApplication
import pickle

glutInit()
app = QApplication(sys.argv)
from pele.gui.neb_explorer import NEBExplorer
system = OXDNASystem()

wnd = NEBExplorer(None, system, app)
wnd.show()
path = pickle.load(open("interpolate.pickle"))
wnd.new_neb(path[0], path[-1], path=path, run=False)

sys.exit(app.exec_())