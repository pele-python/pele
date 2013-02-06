import sys
from cStringIO import StringIO

from PyQt4 import QtGui, QtCore, Qt
from PyQt4.QtGui import QDialog, QApplication, QListWidgetItem


from pygmin.gui.connect_run_ui import Ui_MainWindow as UI
from pygmin.storage import Database

class OutLog(object):
    def __init__(self, edit, out=None, color=None):
        """(edit, out=None, color=None) -> can write stdout, stderr to a
        QTextEdit.
        edit = QTextEdit
        out = alternate stream ( can be the original sys.stdout )
        color = alternate color (i.e. color stderr a different color)
        
        from http://www.riverbankcomputing.com/pipermail/pyqt/2009-February/022025.html
        """
        self.edit = edit
        self.out = None
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


class ConnectViewer(QtGui.QMainWindow):#QtGui.QMainWindow):
    def __init__(self, parent=None, system=None, app=None, decrunner=None):
        QtGui.QMainWindow.__init__(self)
        self.decrunner = decrunner
    
        self.ui = UI()
        self.ui.setupUi(self)
        self.ui.centralwidget.hide()
        self.textEdit = self.ui.textEdit
        self.textEdit.setReadOnly(True)
        self.textEdit.insertPlainText("first line\n")
        self.textEdit.insertPlainText("second line\n")
        
        self.textEdit_writer = OutLog(self.textEdit, sys.stdout)


#
# only testing stuff below here
#


def start():
    decrunner.start()
    print >> sys.stderr, "started decrunner"
#    wnd.finish()


if __name__ == "__main__":
    import sys
    import pylab as pl
    from pygmin.gui.double_ended_connect_runner import DECRunner

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
    
    
    wnd = ConnectViewer(app=app, system=system)
    decrunner = DECRunner(system, db, min1, min2, outstream=wnd.textEdit_writer)
    wnd.show()
    from PyQt4.QtCore import QTimer
    QTimer.singleShot(10, start)
    sys.exit(app.exec_()) 