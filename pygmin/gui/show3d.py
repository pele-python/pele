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
from pygmin.utils.events import Signal
from PyQt4.QtCore import pyqtSlot
#from pygmin.gui.ui.show3d_with_slider_ui import Ui_Form

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
    """
    The main OpenGL viewer for pygmin
    
    This viewer can display a single structure, two structures
    overlaid, or a pathway of structures (with a slider) 
    """
    def __init__(self, *args, **kwargs):
        super(Show3DWithSlider, self).__init__(*args, **kwargs)
        
#        self.ui = Ui_Form()
#        self.ui.setupUi(self)
         
        self.label = QtGui.QLabel(parent=self)
        self.label.setObjectName(_fromUtf8("label"))
        self.label.setText("")

        
        self.oglwgt = Show3D(parent=self)
        
        self.slider = QtGui.QSlider(parent=self)
        self.slider.setOrientation(QtCore.Qt.Horizontal)
        self.slider.setObjectName(_fromUtf8("slider"))
        QtCore.QObject.connect(self.slider, QtCore.SIGNAL(_fromUtf8("valueChanged(int)")), self.on_slider_valueChanged)

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.label)
        vbox.addWidget(self.oglwgt)
        vbox.addWidget(self.slider)
        self.setLayout(vbox)
        
        self.on_frame_updated = Signal()

        self.animate = False
        self._animate_dir = 1
    
#    def __getattr__(self, name):
#        return getattr(self.oglwgt, name)
    def setSystem(self, system):
        """
        set the system which has information about how to draw the structures.
        
        it should have a function:
        
            system.draw(coords)
        """
        self.oglwgt.setSystem(system)
        
    def setCoords(self, coords, index=1):
        """
        display a structure (not a path) in the ogl viewer
        
        Parameters
        ----------
        coords : 1d array
            the structure to display
        index : 1 or 2
            the viewer can display two structures overlayed.  This specifies
            whether this structure should be in slot 1 or slot 2
        """
        if index not in (1, 2):
            raise ValueError("index must be either 1 or 2")
        self.animate = False
        self.messages = None
        self.coordspath = None
        self.slider.hide()
        self.oglwgt.setCoords(coords, index=index)
        self.label.setText("")

#    def getCoords(self, index=1):
#        if index not in (1, 2):
#            raise ValueError("index must be either 1 or 2")
#        return self.oglwgt.coords[index]
#    
#    def setMinimum(self, m, index=1):
#        if index not in (1, 2):
#            raise ValueError("index must be either 1 or 2")
#        self.oglwgt.minima[index] = m
#
#    def getMinimum(self, index=1):
#        if index not in (1, 2):
#            raise ValueError("index must be either 1 or 2")
#        return self.oglwgt.minima[index]
    
    def setCoordsPath(self, coordspath, frame=None, labels=None):
        """
        show a path in the viewer
        
        Parameters
        ----------
        coordspath : 2d numpy array
            coordspath[i,j] is the coordinates of the j'th degree of freedom of the i'th structure
        frame : int
            frame number to show initially
        labels : list of strings
            labels for the structures that will be shown above the ogl viewer
        """
        self.oglwgt.setCoords(None, index=2)
        if(frame is None):
            frame = self.slider.value()
            
        if labels is None:
            self.label.hide()
        else:
            self.label.show()
        
        if frame < 0:
            frame = coordspath.shape[0]-1
        else:
            frame = min(frame, coordspath.shape[0]-1)
        
        self.coordspath = coordspath
        self.messages = labels
        self.slider.show()
        self.slider.setRange(0, coordspath.shape[0]-1)
        self.showFrame(frame)
    
    @pyqtSlot(int)
    def on_slider_valueChanged(self, i):
        print i
        return self._showFrame(i)

    def _showFrame(self, i):
        self.oglwgt.setCoords(self.coordspath[i,:], index=1)
        if self.messages is not None:
            self.label.setText(self.messages[i])
        self.on_frame_updated(i, sender=self)

    def showFrame(self, i):
        if i == -1:
            i = self.coordspath.shape[0] - 1
        self.slider.setValue(i)
        self._showFrame(i)

    def get_slider_index(self):
        return self.slider.value()

#    def setCoordsSingle(self):

#    def on_myslider_sliderMoved(self, index):
#        print "slider moved", index
##        self.oglwgt.setCoords(self.neb.coords[0,:], index=index)


    def start_animation(self):
        self.animate=True
        self._animate_dir = 1
        QtCore.QTimer.singleShot(0., self._next_frame)
    
    def stop_animation(self):
        self.animate = False
    
    def _next_frame(self):
        if not self.animate: return
        cur = self.slider.value()
        if cur == self.slider.maximum():
            self._animate_dir = -1
        elif cur ==  self.slider.minimum():
            self._animate_dir = 1
        cur += self._animate_dir
#        self.slider.setValue(cur)
        self.showFrame(cur)
        
        if self.animate:
            frames_per_second = 10.
            QtCore.QTimer.singleShot(1000. / frames_per_second, self._next_frame)

