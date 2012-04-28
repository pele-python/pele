import numpy as np
import rotations as rot

def vec2aa( v2, v1 = np.array( [0.,0.,1.]) ):
    """ 
    return the angle-axis rotation which transforms
    v1 into v2
    """
    #find perpendicular
    vperp = -np.cross(v1,v2)
    aa = vperp / np.linalg.norm(vperp) * np.dot(v1,v2)
    return aa
    

class Site:
    """
    this class defines a site in a rigid body
    """
    def __init__(self, type, position):
        self.type = type
        assert( len(position) == 3 )
        self.position = np.array(position)
        #self.aa = vec2aa( position )
        #self.rotation_matrix = rot.aa2mx(self.aa)
        


class Molecule:
    """
    this class defines a rigid body
    """
    def __init__(self):
        self.sitelist=[]
        self.nsites = 0
        pass
    
    def insert_site(self, type, position):
        self.sitelist.append( Site(type, position) )
        self.nsites += 1
    
    def getxyz(self, aa=np.array([0.,0.,1e-6])):
        nsites = len(self.sitelist)
        xyz = np.zeros(nsites*3, np.float64)
        mx = rot.aa2mx(aa)
        for i,site in enumerate(self.sitelist):
            xyz[i*3:i*3+3] = np.dot(mx, site.position)
        return xyz

def setupLWOTP():
    from numpy import sin, cos, pi
    pos1 = [0.0, -2./3 * np.sin( 7.*pi/24.), 0.0]
    pos2 = [cos( 7.*pi/24.),  1./3. * sin( 7.* pi/24.), 0.0]
    pos3 = [-cos( 7.* pi/24.),  1./3. * sin( 7.*pi/24), 0.0]
    otp= Molecule()
    otp.insert_site(0, pos1 )
    otp.insert_site(0, pos2 )
    otp.insert_site(0, pos3 )
    return otp
                
            
def test_molecule():
    from numpy import sin, cos, pi
    otp = setupLWOTP()
        
    xyz = otp.getxyz()
    from printing.print_atoms_xyz import printAtomsXYZ as printxyz
    import sys
    #with open("out.xyz", "w") as fout:
    printxyz(sys.stdout, xyz)
    
    
    
if __name__ == "__main__":
    test_molecule()