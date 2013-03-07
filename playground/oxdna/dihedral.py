import numpy as np
from numpy.linalg import norm
import numpy as np

def dihedral_angle(r):
        """ r1,r2,r3,r4 are numpy array containing x,y,z coordinates

        return angle is between 0 and 360 degrees
        """
        r1 = r[0]
        r2 = r[1]
        r3 = r[2]
        r4 = r[3]

        # normal to plane 1-2-3
        normal1 = np.cross( r2 - r1 , r3 - r2 )
        normal1 = normal1 / norm(normal1)

        # normal to plane 2-3-4
        normal2 = np.cross( r3 - r2 , r4 - r3 )
        normal2 = normal2 / norm(normal2)

        # cos ( angle between normals )
        costheta = np.dot( normal1, np.transpose(normal2) )

        # check if cross product of normals is parallel or antiparallel
        # to vector r3-r2 connecting two planes
        anchor =  ( r3 - r2) / norm(r3 -r2)
        nnormal = np.cross(normal1, normal2)

        cosNormal = np.dot( anchor, np.transpose(nnormal))

        if costheta > 1.0:
            costheta = 1.0
        if costheta < -1.0:
            costheta = -1.0
        rad = np.arccos(costheta)

        if cosNormal < 0:
            rad = 2*np.pi - np.arccos(costheta)

        return rad

def dihedral_gradient(r):
        """ r1,r2,r3,r4 are numpy array containing x,y,z coordinates
        returns derivative of angle wrt particle position)"""
        r1 = r[0]
        r2 = r[1]
        r3 = r[2]
        r4 = r[3]
        g = np.zeros([4,3])

        r12 = r1 - r2
        r32 = r3 - r2
        r34 = r3 - r4

        b12 = np.linalg.norm(r12)
        b32 = np.linalg.norm(r32)
        b34 = np.linalg.norm(r34)

        b32s = b32**2
        
        normal_a = np.cross( r12 , r32 )
        normal_b = np.cross( r32 , r34 )

        deltaX_deltar1 = b32 / ( np.dot( normal_a,normal_a )) * normal_a

        deltaX_deltar4 = - b32 / ( np.dot( normal_b, normal_b )) * normal_b

        deltaX_deltar2 = ((np.dot(r12,r32) / b32s) - 1) \
            * deltaX_deltar1 - ((np.dot(r34,r32) / b32s)) \
            * deltaX_deltar4

        deltaX_deltar3 = ((np.dot(r34,r32) / b32s) - 1) \
            * deltaX_deltar4 - ((np.dot(r12,r32) / b32s)) \
            * deltaX_deltar1

        g[0,:] =  deltaX_deltar1
        g[1,:] =  deltaX_deltar2
        g[2,:] =  deltaX_deltar3
        g[3,:] =  deltaX_deltar4
        
        return g


