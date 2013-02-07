import sys
import numpy as np

from PyQt4 import QtGui, QtCore, Qt
from PyQt4.QtGui import QDialog, QApplication, QListWidgetItem


from pygmin.gui.connect_run_ui import Ui_MainWindow as UI
from pygmin.utils.events import Signal
from pygmin.gui.double_ended_connect_runner import DECRunner


class OutLog(object):
    """recieve text through self.write and write it to a QtextEdit object
    
    (edit, out=None, color=None) -> can write stdout, stderr to a
    QTextEdit.
    
    Parameters
    ----------
    edit : QTextEdit
    out : something like a filestream (has out.write(mystring))
        alternate stream ( can be the original sys.stdout )
    color : 
        alternate color (i.e. color stderr a different color)
    
    from http://www.riverbankcomputing.com/pipermail/pyqt/2009-February/022025.html
    """
    def __init__(self, edit, out=None, color=None):
        self.edit = edit
        self.out = out
        self.color = color

    def write(self, m):
        if self.color:
            tc = self.edit.textColor()
            self.edit.setTextColor(self.color)

        self.edit.moveCursor(QtGui.QTextCursor.End)
        self.edit.insertPlainText( m )

        if self.color:
            self.edit.setTextColor(tc)

        if self.out:
            self.out.write(m)
    
    def flush(self):
        pass


class ConnectViewer(QtGui.QMainWindow):
    """
    external viewer for connect runs
    
    This viewer will use DECRunner to run a connect job in parallel.
    The log messages from that run will be redirected into a GUI text viewer.
    When the run is finished, the completed path will be shown in an OGL viewer
    
    Parameters
    ----------
    system : 
    database :
    min1, min2 : Minimum objects
        the minima to try to connect
    parent : 
        the parent window
    app : 
        the application
    
    See Also
    --------
    DECRunner
    
    """
    def __init__(self, system, database, min1, min2, parent=None, app=None):
        QtGui.QMainWindow.__init__(self, parent=parent)    
        self.ui = UI()
        self.ui.setupUi(self)
        self.ui.centralwidget.hide()
        
        self.ogl = self.ui.ogl
        self.ogl.setSystem(system)

        
        self.textEdit = self.ui.textEdit      
        self.textEdit.setReadOnly(True)
        self.textEdit_writer = OutLog(self.textEdit)
        
        self.decrunner = DECRunner(system, database, min1, min2, outstream=self.textEdit_writer)
        self.decrunner.on_finished.connect(self.on_finished)

    def start(self):
        self.decrunner.start()

    def on_finished(self):
        print "success", self.decrunner.success
#        print "success", self.decrunner.smoothed_path
        if self.decrunner.success:
            self.smoothed_path = np.array(self.decrunner.smoothed_path)
#            print self.smoothed_path.shape
            self.ogl.setCoordsPath(self.smoothed_path)

        

#
# only testing stuff below here
#


def start():
    wnd.start()
    print >> sys.stderr, "started decrunner"
#    wnd.finish()


if __name__ == "__main__":
    from OpenGL.GLUT import glutInit
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
    db = system.create_database()
    min1 = db.addMinimum(e1, x1)
    min2 = db.addMinimum(e2, x2)
    
    
    wnd = ConnectViewer(system, db, min1, min2, app=app)
#    decrunner = DECRunner(system, db, min1, min2, outstream=wnd.textEdit_writer)
    glutInit()
    wnd.show()
    from PyQt4.QtCore import QTimer
    QTimer.singleShot(10, start)
    sys.exit(app.exec_()) 