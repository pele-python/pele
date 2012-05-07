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
from PyQt4.QtOpenGL import *
import numpy as np

class Show3D(QGLWidget):
    '''
    Widget for drawing two spirals.
    '''
    
    def __init__(self, parent):
        QGLWidget.__init__(self, parent)
        self.setMinimumSize(200, 200)
        glutInit()#sys.argv)
        self.coords = None
        self.last_mouse_pos=QtCore.QPointF(0.,0.)
        self.rotation=np.array([0.,0.])

    def setCoords(self, coords):
        self.coords = coords.reshape(coords.size/3,3)
        
    def paintGL(self):
        '''
        Drawing routine
        '''
         
    	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
    	if(self.coords == None):
            return
       
        color = [0.65,0.0,0.0,1.]
        amb = [0.,0.0,0.0,1.]
        spec=[1.0,1.,1.,1]
        shine=80.
        glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,color)
        glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,amb)
        glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,shine)
        glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,spec)
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glRotate(self.rotation[0], 1., 0.,0.)
        glRotate(self.rotation[1], 0., 1.,0.)
        com=np.mean(self.coords, axis=0)      
            
        for xx in self.coords:
            x=xx-com
            glPushMatrix()            
            glTranslate(x[0],x[1],x[2])
    	    glutSolidSphere(0.5,30,30)
    	    glPopMatrix()
        glPopMatrix()	
        glFlush()
    	#glutSwapBuffers()

    def mouseMoveEvent(self, event):
        self.last_mouse_pos
        delta = event.posF() - self.last_mouse_pos
        self.last_mouse_pos=event.posF()
        self.rotation+=np.array([delta.y(), delta.x()])
        self.repaint()
        
    def resizeGL(self, w, h):
        '''
        Resize the GL window 
        '''
        print "resizegl"
        glViewport(0, 0, w, h)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(40.0, 1.0, 1.0, 40.0)
    
    def initializeGL(self):
        '''
        Initialize GL
        '''
        
        glClearColor(1.,1.,1.,1.)
    	glShadeModel(GL_SMOOTH)
    	glEnable(GL_CULL_FACE)
    	glEnable(GL_DEPTH_TEST)
    	glEnable(GL_LIGHTING)
    	lightZeroPosition = [-80.,50.,150.,1.]
    	lightZeroColor = [1.0, 1.,1.,1.0] #green tinged
    	glLightfv(GL_LIGHT0, GL_POSITION, lightZeroPosition)
    	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor)
        glLightfv(GL_LIGHT0, GL_AMBIENT, [0.5,0.5,0.5,1.0])
        
        
        glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 1.000)
    	#glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.01)
    	glEnable(GL_LIGHT0)
    	
        glLightfv(GL_LIGHT1, GL_POSITION, [10.,20.,-2.,1.])
        glLightfv(GL_LIGHT1, GL_DIFFUSE, [1.,1.,1.,1.0])
        glLightfv(GL_LIGHT1, GL_SPECULAR, [1.0,1.0,1.0,1.0])
        glLightfv(GL_LIGHT1, GL_AMBIENT, [0.,0.,0.,1.0])
        glLightf(GL_LIGHT1, GL_CONSTANT_ATTENUATION, 1)
        #glLightf(GL_LIGHT1, GL_LINEAR_ATTENUATION, 0.01)
        glEnable(GL_LIGHT1)
        
        glShadeModel(GL_SMOOTH);
        glMatrixMode(GL_PROJECTION)
    	gluPerspective(40.,1.,1.,40.)
    	glMatrixMode(GL_MODELVIEW)
    	gluLookAt(0,0,10,
    		  0,0,0,
    		  0,1,0)
    	glPushMatrix()
    	
	
