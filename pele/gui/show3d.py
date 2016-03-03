"""
Created on Jul 7, 2009

@author: Stou Sandalski (stou@icapsid.net)
@license:  Public Domain
"""
import numpy as np

from OpenGL import GL
from OpenGL.GL import glMaterialfv, glEnable, glLightfv
from OpenGL import GLU
from PyQt4 import QtGui, QtCore
from PyQt4.Qt import Qt
from PyQt4.QtOpenGL import QGLWidget

import pele.utils.rotations as rot

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s
try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)



class Show3D(QGLWidget):
    """
    widget for displaying OpenGL.  This exists only for legacy support.  use Show3DWithSlider instead
    """
    
    def __init__(self, parent):
        QGLWidget.__init__(self, parent)
        self.setMinimumSize(200, 200)
#        glutInit()#sys.argv)
        self.coords = {1:None, 2:None}
        self.minima = {1:None, 2:None}
        self.last_mouse_pos = QtCore.QPointF(0., 0.)
        self.rotation = rot.aa2mx(np.array([0.,0.,0.])) # np.array([0., 0.])
        self.zoom = 1.0
        self._fatal_error = False # don't try to plot if it won't work
        
    def setCoords(self, coords, index=1):
        self.coords[index] = coords
        self.repaint()

    def setMinimum(self, minimum, index=1):
        self.minima[index] = minimum
        self.repaint()
        
    def setSystem(self, system):
        self.system = system
        
    def paintGL(self):
        """
        Drawing routine
        """
        if self._fatal_error: return
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        GL.glEnable(GL.GL_NORMALIZE)
        amb = [0., 0.0, 0.0, 1.]
        spec = [1.0, 1., 1., 1]
        shine = 80.
        glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_AMBIENT, amb)
        glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_SHININESS, shine)
        glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_SPECULAR, spec)
        GL.glMatrixMode(GL.GL_MODELVIEW)
        GL.glPushMatrix()
        GL.glScalef(self.zoom, self.zoom, self.zoom)
        #glRotate(self.rotation[0], 1., 0., 0.)
        #glRotate(self.rotation[1], 0., 1., 0.)
        mx=np.zeros([4,4])
        mx[0:3,0:3] = self.rotation
        mx[3,3]=1.
        GL.glMultMatrixf(mx)
        
        for index, coords in self.coords.items():
            if coords is None:
                continue
            if index == 1:
                color = [0.65, 0.0, 0.0, 1.]
            if index == 2:
                color = [0.00, 0.65, 0., 1.]
            GL.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE, color)
            try:
                self.system.draw(coords, index)
            except NotImplementedError:
                self._fatal_error = True
                raise
        
        GL.glPopMatrix()    
        GL.glFlush()
        #glutSwapBuffers() @IndentOk

    def mousePressEvent(self, event):
        self.last_mouse_pos = event.posF()
    
    def mouseMoveEvent(self, event):
        delta = (event.posF() - self.last_mouse_pos)*0.01
        self.last_mouse_pos = event.posF()
        if event.buttons() == Qt.LeftButton:
            drot = rot.aa2mx(-np.array([delta.y(), delta.x(), 0.]))
            self.rotation = np.dot(self.rotation, drot)
        elif event.buttons() == Qt.RightButton:
            drot = rot.aa2mx(np.array([0., 0., delta.x()]))
            self.rotation = np.dot(self.rotation, drot)
            self.zoom *= 1.0 - 0.2*delta.y()
        self.repaint()
        
    def resizeGL(self, w, h):
        """
        Resize the GL window
        """
        GL.glViewport(0, 0, w, h)
        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glLoadIdentity()
        GLU.gluPerspective(40.0, float(w)/float(h), 1.0, 40.0)
    
    def initializeGL(self):
        """
        Initialize GL
        """
        
        GL.glClearColor(1., 1., 1., 1.)
        GL.glShadeModel(GL.GL_SMOOTH)
        glEnable(GL.GL_CULL_FACE)
        glEnable(GL.GL_DEPTH_TEST)
        glEnable(GL.GL_LIGHTING)
        lightZeroPosition = [-800., 500., 1500., 1.]
        lightZeroColor = [1.0, 1., 1., 1.0] #green tinged
        glLightfv(GL.GL_LIGHT0, GL.GL_POSITION, lightZeroPosition)
        glLightfv(GL.GL_LIGHT0, GL.GL_DIFFUSE, lightZeroColor)
        glLightfv(GL.GL_LIGHT0, GL.GL_AMBIENT, [0.5, 0.5, 0.5, 1.0])
        
        
        GL.glLightf(GL.GL_LIGHT0, GL.GL_CONSTANT_ATTENUATION, 1.000)
        #glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.01)
        glEnable(GL.GL_LIGHT0)
        
        glLightfv(GL.GL_LIGHT1, GL.GL_POSITION, [100., 200., -20., 1.])
        glLightfv(GL.GL_LIGHT1, GL.GL_DIFFUSE, [1., 1., 1., 1.0])
        glLightfv(GL.GL_LIGHT1, GL.GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
        glLightfv(GL.GL_LIGHT1, GL.GL_AMBIENT, [0., 0., 0., 1.0])
        GL.glLightf(GL.GL_LIGHT1, GL.GL_CONSTANT_ATTENUATION, 1)
        #glLightf(GL_LIGHT1, GL_LINEAR_ATTENUATION, 0.01)
        glEnable(GL.GL_LIGHT1)
        
        GL.glShadeModel(GL.GL_SMOOTH)
        GL.glMatrixMode(GL.GL_PROJECTION)
        GLU.gluPerspective(40., 1., 1., 40.)
        GL.glMatrixMode(GL.GL_MODELVIEW)
        GLU.gluLookAt(0, 0, 10,
              0, 0, 0,
              0, 1, 0)
        GL.glPushMatrix()
        
    def sizeHint(self):
        w, h = 500,500 #self.get_width_height()
        return QtCore.QSize(w, h)

    def minimumSizeHint(self):
        return QtCore.QSize(10, 10)

