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
import pygmin.utils.rotations as rot

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s


class Show3D(QGLWidget):
    '''
    Widget for drawing two spirals.
    '''
    
    def __init__(self, parent):
        QGLWidget.__init__(self, parent)
        self.setMinimumSize(200, 200)
#        glutInit()#sys.argv)
        self.coords = {}
        self.minima = {}
        self.last_mouse_pos = QtCore.QPointF(0., 0.)
        self.rotation = rot.aa2mx(np.array([0.,0.,0.])) # np.array([0., 0.])
        self.zoom = 1.0
        
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
            self.system.draw(coords, index)
        
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
        print "resizegl"
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

class Show3DWithSlider(QWidget):
    def __init__(self, *args, **kwargs):
        super(Show3DWithSlider, self).__init__(*args, **kwargs)
        
        self.oglwgt = Show3D(parent=self)
        
        self.slider = QtGui.QSlider(parent=self)
        self.slider.setOrientation(QtCore.Qt.Horizontal)
        self.slider.setObjectName(_fromUtf8("myslider"))
        QtCore.QObject.connect(self.slider, QtCore.SIGNAL(_fromUtf8("sliderMoved(int)")), self._showFrame)

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.oglwgt)
        vbox.addWidget(self.slider)
        self.setLayout(vbox)

    
#    def __getattr__(self, name):
#        return getattr(self.oglwgt, name)
    def setSystem(self, system):
        self.oglwgt.setSystem(system)
    
    def setCoordsPath(self, coordspath, frame=0):
        self.coordspath = coordspath
        self.slider.setRange(0, coordspath.shape[0]-1)
        self.showFrame(frame)

    def _showFrame(self, i):
        self.oglwgt.setCoords(self.coordspath[i,:], index=1)
        

    def showFrame(self, i):
        self.slider.setValue(i)
        self._showFrame(i)

#    def setCoordsSingle(self):

#    def on_myslider_sliderMoved(self, index):
#        print "slider moved", index
##        self.oglwgt.setCoords(self.neb.coords[0,:], index=index)

