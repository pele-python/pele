from OpenGL.GL import *
from OpenGL.GLU import *
from PyQt4 import QtGui
from PyQt4.QtOpenGL import *
from PyQt4.Qt import Qt
from PyQt4 import QtCore

from . import pymol2




class PymolQtWidget(QGLWidget):
     _buttonMap = {Qt.LeftButton:0,
                 Qt.MidButton:1,
                 Qt.RightButton:2}


     def __init__(self, parent, enableUi,File=""):
         f = QGLFormat()
         f.setStencil(True)
         f.setRgba(True)
         f.setDepth(True)
         f.setDoubleBuffer(True)
         QGLWidget.__init__(self, f, parent=parent)
         self.setMinimumSize(500, 500)
         self._enableUi=enableUi
         self.pymol = pymol2.PyMOL()# _pymolPool.getInstance()
         self.pymol.start()
         self.cmd = self.pymol.cmd
         # self.toPymolName = self.pymol.toPymolName ### Attribute Error
         self._pymolProcess()
         if not self._enableUi:
             self.pymol.cmd.set("internal_gui",0)
             self.pymol.cmd.set("internal_feedback",0)
             self.pymol.cmd.button("double_left","None","None")
             self.pymol.cmd.button("single_right","None","None")

         self.pymol.cmd.load(File)
         self.pymol.reshape(self.width(),self.height())
         self._timer = QtCore.QTimer()
         self._timer.setSingleShot(True)
         self._timer.timeout.connect(self._pymolProcess)
         self.resizeGL(self.width(),self.height())
         #globalSettings.settingsChanged.connect(self._updateGlobalSettings)
         self._updateGlobalSettings()

     def __del__(self):
         pass

     def _updateGlobalSettings(self):
         #for k,v in globalSettings.settings.iteritems():
         #    self.pymol.cmd.set(k, v)
         #self.update()
         return

     def redoSizing(self):
         self.resizeGL(self.width(), self.height())

     def paintGL(self):
         glViewport(0,0,self.width(), self.height())
         bottom = self.mapToGlobal(QtCore.QPoint(0,self.height())).y()
         #self.pymol.cmd.set("_stencil_parity", bottom & 0x1)
         self._doIdle()
         self.pymol.draw()

     def mouseMoveEvent(self, ev):
         self.pymol.drag(ev.x(), self.height()-ev.y(),0)
         self._pymolProcess()

     def mousePressEvent(self, ev):
         if not self._enableUi:
             self.pymol.cmd.button("double_left","None","None")
             self.pymol.cmd.button("single_right","None","None")
         self.pymol.button(self._buttonMap[ev.button()], 0, ev.x(),
self.height()-ev.y(),0)
         self._pymolProcess()

     def mouseReleaseEvent(self, ev):
         self.pymol.button(self._buttonMap[ev.button()], 1, ev.x(),
self.height()-ev.y(),0)
         self._pymolProcess()
         self._timer.start(0)

     def resizeGL(self, w, h):
         self.pymol.reshape(w,h, True)
         self._pymolProcess()

     def initializeGL(self):
         pass

     def _pymolProcess(self):
         self._doIdle()
         self.update()

     def _doIdle(self):
         if self.pymol.idle():
             self._timer.start(0)



# You don't need anything below this
class PyMolWidgetDemo(QtGui.QMainWindow):
     def __init__(self):
         QtGui.QMainWindow.__init__(self)
         widget = PymolQtWidget(self,True,"D2.xyz")
         self.setCentralWidget(widget)

if __name__ == '__main__':
     app = QtGui.QApplication(['PyMol Widget Demo'])
     window = PyMolWidgetDemo()
     window.show()
     app.exec_()

