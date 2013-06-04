
import numpy as np

__all__ = ["Measure"]

class Measure:
    """
    Measure length, angle and torsion angle given Cartesian points in 3-D
    
    """ 
    def __init__(self):
        self.dummy = 0 
    
    def norm(self, r1):
        """ norm of a vector  """        
        return np.sqrt( np.sum( np.square( r1 )))

    def bond(self, r1, r2):
        """ r1, r2 are numpy array with x,y,z coordinates """        
        return np.sqrt( np.sum( np.square( r1-r2 )))
        
    def angle(self, r1, r2, r3):
        """ r1,r2,r3 are numpy array with x,y,z coordinates 
        
        return angle is between 0 and 180 degrees  
        """
        
        r21  =  r1 - r2
        r23  =  r3 - r2 
        
        b12 = self.bond(r1, r2) 
        b23 = self.bond(r2, r3)
                
        costheta = np.divide( np.dot(np.transpose(r21),r23) , np.dot(b12,b23))                      

        rad = np.arccos(costheta)   
        deg = rad*180/np.pi  
            
        return rad, deg 

    def torsion(self, r1, r2, r3, r4):
        """ r1,r2,r3,r4 are numpy array containing x,y,z coordinates 
        
        return angle is between 0 and 360 degrees  
        """
        
        # normal to plane 1-2-3
        normal1 = np.cross( r2 - r1 , r3 - r2 )  
        normal1 = normal1 / self.norm(normal1) 
        
        # normal to plane 2-3-4
        normal2 = np.cross( r3 - r2 , r4 - r3 )  
        normal2 = normal2 / self.norm(normal2) 
         
        # cos ( angle between normals ) 
        costheta = np.dot( normal1, np.transpose(normal2) )  
        
        # check if cross product of normals is parallel or antiparallel 
        # to vector r3-r2 connecting two planes 
        anchor =  ( r3 - r2) / self.norm(r3 -r2)  
        nnormal = np.cross(normal1, normal2) 
        
        cosNormal = np.dot( anchor, np.transpose(nnormal))   
        
        rad = np.arccos(costheta)  
        
        if cosNormal < 0: 
            rad = 2*np.pi - np.arccos(costheta)   

        deg = rad*180/np.pi           ;

        return rad, deg 

if __name__ == "__main__":
    
    a = np.array([1,2,3])
    b = np.array([0,2,3])

    bat = Measure() 
         
    print 'bond length (should be 1) = '         
    print bat.bond(a,b)
    
    print 'angle (should be 90 deg) = '         
    print bat.angle( np.array([1, 0 , 0 ]) , np.array([0, 0 ,0 ]) , np.array([0, -1, 0]))

    r1 = np.array([ 1, 0, 0 ]) #  on x axis
    r2 = np.array([ 0, 0, 0 ]) #  origin 
    r3 = np.array([ 0, 1, 0 ]) #  on y axis 
    r4 = np.array([ 0, 1, 1 ]) #  in y-z plane 
    print 'torsion angle (should be 270 deg) = '         
    print bat.torsion( r1, r2, r3, r4)

    print 'all done'
    
    
    