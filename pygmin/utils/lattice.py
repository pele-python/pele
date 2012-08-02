'''Helper functions to handle crystal lattice'''
import numpy as np
import collections
from math import sqrt, cos, sin, acos

def lowerTriangular(coords):
    ''' convert 6 lattice degrees of freedom to lower triangular lattice matrix
    
        :param coords: 6 lattice coordinates
        :returns: lattice matrix    
    '''
    
    return np.array(        
        [[coords[0],    0.,        0.],
        [ coords[1], coords[3],    0.],
        [ coords[2], coords[4], coords[5]]])

def volume(coords):
    ''' calcualte box volume
        :param coords: 6 lattice coordinates
        :returns: box volume
    '''
    
    ml = lowerTriangular(coords)
    return np.dot(ml[:,0],np.cross(ml[:,1], ml[:,2]))

class Cell(object):
    ''' class to handle crystal cells 
    
        The internal representation of the lattice is in the gruber notation.
    '''
    
    _a = 0.
    _b = 0.
    _c = 0.
    _d = 0.
    _e = 0.
    _f = 0.
    
    def set_lattice_matrix(self, lattice):
        ''' set cell parameters from lattice matrix
        
            :param lattice: 3x3 lattice matrix            
        '''
        la = lattice[:,0]
        lb = lattice[:,1]
        lc = lattice[:,2]

        self._a = np.dot(la,la)
        self._b = np.dot(lb,lb)
        self._c = np.dot(lc,lc)
        self._d = 2.0*np.dot(lb,lc)
        self._e = 2.0*np.dot(la,lc)
        self._f = 2.0*np.dot(la,lb)
    
    
    def set_with_angles(self, a, b, c, alpha, beta, gamma):
        ''' set cell parameters from angle.
            
            :param a: length of box vector a
            :type a: float
            :param b: length of box vector b
            :type b: float
            :param c: length of box vector c
            :type c: float
            :param alpha: angle between b and c in rad
            :type alpha: float
            :param beta: angle between a and c in rad
            :type beta: float
            :param gamma: angle between a and b in rad
            :type gamma: float
        '''
        self._a = a*a
        self._b = b*b
        self._c = c*c
        self._d = 2.0*cos(alpha)*b*c
        self._e = 2.0*cos(beta)*a*c
        self._f = 2.0*cos(gamma)*a*b
        
    def get_with_angles(self):
        ''' get the cell parameters with box length / angles
            :returns: named tuple a,b,c,alpha,beta,gamma
        '''
        a,b,c = sqrt(self._a), sqrt(self._b), sqrt(self._c)
        return collections.namedtuple("CellAngle", "a, b, c, alpha, beta, gamma")(
           a, b, c,
           acos(0.5 * self._d / b / c),
           acos(0.5 * self._e_ / a / c),
           acos(0.5 * self._f_ / a / b))
    
    def get_lower_triangular(self):
        ''' calculate the lower triangular lattice matrix
            :returns: lower triangular matrix
        '''
        p = self.get_with_angles()
        m=np.zeros(3,3)
        m[2,2]=p.c
        m[1,1]=p.b*sin(p.alpha)
        m[2,1]=p.b*cos(p.alpha)
        m[0,0]=p.a/sin(p.alpha)*sqrt(1. - cos(p.alpha)**2 - cos(p.beta)**2 \
                 - cos(p.gamma)**2  + 2*cos(p.alpha)*cos(p.beta)*cos(p.gamma))
        m[1,0]=p.a*(cos(p.gamma) - cos(p.beta)*cos(p.alpha))/sin(p.alpha)
        m[2,0]=p.a*cos(p.beta)
        return m
    
    def get_upper_triangular(self):
        ''' calculate the upper triangular lattice matrix
            :returns: upper triangular matrix
        '''
        p = self.get_with_angles()
        m=np.zeros(3,3)
        
        v=sqrt(1-cos(p.alpha)*cos(p.alpha)-cos(p.beta)*cos(p.beta)\
               -cos(p.gamma)*cos(p.gamma)+2*cos(p.alpha)*cos(p.beta)*cos(p.gamma))

        #aa=sin(p.alpha)/p.a/v
        #bb=sin(p.beta )/p.b/v
        cc=sin(p.gamma)/p.c/v

        alphaa=acos( (cos(p.beta )*cos(p.gamma)-cos(p.alpha))/sin(p.beta )/sin(p.gamma) )
        #betaa =acos( (cos(p.alpha)*cos(p.gamma)-cos(p.beta ))/sin(p.alpha)/sin(p.gamma) )
        #gammaa=acos( (cos(p.alpha)*cos(p.beta )-cos(p.gamma))/sin(p.alpha)/sin(p.beta ) )

        m[0,0]=p.a
        m[0,1]=p.b*cos(p.gamma)
        m[0,2]=p.c*cos(p.beta)
        
        m[1,0]=0
        m[1,1]=p.b*sin(p.gamma)
        m[1,2]=-p.c*sin(p.beta)*cos(alphaa)
     
        m[2,0]=0
        m[2,1]=0
        m[2,2]=1/cc