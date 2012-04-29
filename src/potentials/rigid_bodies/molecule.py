import numpy as np
import rotations as rot
from potentials.fortran.rmdrvt import rmdrvt as rotMatDeriv

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
        self.abs_position = np.array(position)
        self.drdp = np.zeros([3,3])
        #self.aa = vec2aa( position )
        #self.rotation_matrix = rot.aa2mx(self.aa)



class Molecule:
    """
    this class defines a rigid body molecule
    """
    def __init__(self):
        self.sitelist=[]
        self.nsites = 0

        self.aa = np.zeros(3)
        self.com = np.zeros(3)
        self.rotation_mat = np.zeros([3,3])
        self.drmat = [np.zeros([3,3]) for i in range(3)]
        self.comgrad = np.zeros(3)
        self.aagrad = np.zeros(3)




    def insert_site(self, type, position):
        """
        add a site to the definition of the molecule
        """
        self.sitelist.append( Site(type, position) )
        self.nsites += 1

    def getxyz(self, aa=np.array([0.,0.,1e-6])):
        """return the xyz positions of all sites in the molecule-frame"""
        nsites = len(self.sitelist)
        xyz = np.zeros(nsites*3, np.float64)
        mx = rot.aa2mx(aa)
        for i,site in enumerate(self.sitelist):
            xyz[i*3:i*3+3] = np.dot(mx, site.position)
        return xyz

    def update(self, com, aa, do_derivatives = True):
        """
        Update the position and orientation of the molecule and things dependent on these.
        """
        self.aa[:] = aa
        self.com[:] = com
        self.rotation_mat, self.drmat[0], self.drmat[1], self.drmat[2] = rotMatDeriv(aa, do_derivatives)
        for site in self.sitelist:
            site.abs_position = np.dot( self.rotation_mat, site.position ) + self.com
            for k in range(3):
                site.drdp[k,:] = np.dot(self.drmat[k], site.position)

        #zero interaction dependent things
        self.E = 0.
        self.comgrad[:] = 0.
        self.aagrad[:] = 0.


def setupLWOTP():
    """
    return a molecule of OTP (ortho terphenyl) using the model by Lewis and Wahnstrom
    """
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
