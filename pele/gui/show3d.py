'''
Created on Jul 7, 2009

@author: Stou Sandalski (stou@icapsid.net)
@license:  Public Domain
'''

import math

from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from PyQt4 import QtGui, QtCore
from PyQt4.Qt import Qt, QWidget
from PyQt4.QtOpenGL import *
import numpy as np
import pele.utils.rotations as rot
from pele.utils.events import Signal
from PyQt4.QtCore import pyqtSlot

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
    '''
    widget for displaying OpenGL.  This exists only for legacy support.  use Show3DWithSlider instead
    '''
    
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
        '''
        Drawing routine
        '''
        if self._fatal_error: return
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glEnable(GL_NORMALIZE)
        amb = [0., 0.0, 0.0, 1.]
        spec = [1.0, 1., 1., 1]
        shine = 80.
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, amb)
        glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shine)
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, spec)
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glScalef(self.zoom, self.zoom, self.zoom)
        #glRotate(self.rotation[0], 1., 0., 0.)
        #glRotate(self.rotation[1], 0., 1., 0.)
        mx=np.zeros([4,4])
        mx[0:3,0:3] = self.rotation
        mx[3,3]=1.
        glMultMatrixf(mx)
        
        for index, coords in self.coords.items():
            if coords == None:
                continue
            if index == 1:
                color = [0.65, 0.0, 0.0, 1.]
            if index == 2:
                color = [0.00, 0.65, 0., 1.]
            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color)
            try:
                self.system.draw(coords, index)
            except NotImplementedError:
                self._fatal_error = True
                raise
        
        glPopMatrix()    
        glFlush()
        #glutSwapBuffers() @IndentOk

    def mousePressEvent(self, event):
        self.last_mouse_pos = event.posF()
    
    def mouseMoveEvent(self, event):
        self.last_mouse_pos
        delta = (event.posF() - self.last_mouse_pos)*0.01
        self.last_mouse_pos = event.posF()
        if(event.buttons() == Qt.LeftButton):
            drot = rot.aa2mx(-np.array([delta.y(), delta.x(), 0.]))
            self.rotation = np.dot(self.rotation, drot)
        elif(event.buttons() == Qt.RightButton):
            drot = rot.aa2mx(np.array([0., 0., delta.x()]))
            self.rotation = np.dot(self.rotation, drot)
            self.zoom *= 1.0 - 0.2*delta.y()
        self.repaint()
        
    def resizeGL(self, w, h):
        '''
        Resize the GL window 
        '''
        glViewport(0, 0, w, h)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(40.0, float(w)/float(h), 1.0, 40.0)
    
    def initializeGL(self):
        '''
        Initialize GL
        '''
        
        glClearColor(1., 1., 1., 1.)
        glShadeModel(GL_SMOOTH)
        glEnable(GL_CULL_FACE)
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_LIGHTING)
        lightZeroPosition = [-800., 500., 1500., 1.]
        lightZeroColor = [1.0, 1., 1., 1.0] #green tinged
        glLightfv(GL_LIGHT0, GL_POSITION, lightZeroPosition)
        glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor)
        glLightfv(GL_LIGHT0, GL_AMBIENT, [0.5, 0.5, 0.5, 1.0])
        
        
        glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 1.000)
        #glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.01)
        glEnable(GL_LIGHT0)
        
        glLightfv(GL_LIGHT1, GL_POSITION, [100., 200., -20., 1.])
        glLightfv(GL_LIGHT1, GL_DIFFUSE, [1., 1., 1., 1.0])
        glLightfv(GL_LIGHT1, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
        glLightfv(GL_LIGHT1, GL_AMBIENT, [0., 0., 0., 1.0])
        glLightf(GL_LIGHT1, GL_CONSTANT_ATTENUATION, 1)
        #glLightf(GL_LIGHT1, GL_LINEAR_ATTENUATION, 0.01)
        glEnable(GL_LIGHT1)
        
        glShadeModel(GL_SMOOTH);
        glMatrixMode(GL_PROJECTION)
        gluPerspective(40., 1., 1., 40.)
        glMatrixMode(GL_MODELVIEW)
        gluLookAt(0, 0, 10,
              0, 0, 0,
              0, 1, 0)
        glPushMatrix()
        
    def sizeHint(self):
        w, h = 500,500 #self.get_width_height()
        return QtCore.QSize(w, h)

    def minimumSizeHint(self):
        return QtCore.QSize(10, 10)

