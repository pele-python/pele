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
        self.position = np.array(position) #position in molecule frame
        self.abs_position = np.array(position) #absolute position in space
        self.drdp = np.zeros([3,3]) #derivatives of position w.r.t. aa
        #self.aa = vec2aa( position )
        #self.rotation_matrix = rot.aa2mx(self.aa)
        self.energy = 0.
        self.gradient = np.zeros(3)
        
        self.index = 0 #used to order sites in the system



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
        self.drmat = [np.zeros([3,3]) for i in range(3)] #derivatives of the rotation matrix
        self.comgrad = np.zeros(3)  #gradient of the center of mass
        self.aagrad = np.zeros(3)   #gradient of the angle-axis coords
        
        self.symmetrylist_rot = [] #list of rotational symmetries
        #TODO: implement other types of symmetry



    def insert_site(self, type, position):
        """
        add a site to the definition of the molecule
        """
        self.sitelist.append( Site(type, position) )
        self.nsites += 1
    
    def correctCoM(self):
        """
        perform a few tasks after the sites have been fully defined
        """
        #move the center of mass to the origin
        com = np.zeros(3)
        for site in self.sitelist:
            com += site.position
        com /= self.nsites
        for site in self.sitelist:
            site.position -= com

    def getxyz_rmat(self, rmat, com = np.zeros(3)):
        """return the xyz positions of all sites in the molecule-frame"""
        xyz = np.zeros(self.nsites*3, np.float64)
        for i,site in enumerate(self.sitelist):
            xyz[i*3:i*3+3] = np.dot(rmat, site.position)
        return xyz


    def getxyz(self, com=np.array([0.,0.,0.]), aa=np.array([0.,0.,1e-6]) ):
        """return the xyz positions of all sites in the molecule-frame"""
        nsites = len(self.sitelist)
        xyz = np.zeros(nsites*3, np.float64)
        mx = rot.aa2mx(aa)
        for i,site in enumerate(self.sitelist):
            xyz[i*3:i*3+3] = np.dot(mx, site.position) + com
        return xyz

    def update_coords(self, com, aa, do_derivatives = True):
        """
        Update the position and orientation of the molecule and things dependent on these.
        """
        self.aa[:] = aa
        self.com[:] = com
        self.rotation_mat, self.drmat[0], self.drmat[1], self.drmat[2] = rotMatDeriv(aa, do_derivatives)
        for site in self.sitelist:
            site.abs_position[:] = np.dot( self.rotation_mat, site.position ) + self.com
            for k in range(3):
                site.drdp[k,:] = np.dot(self.drmat[k], site.position)

    def zeroEnergyGrad(self):
        #zero interaction dependent things
        for site in self.sitelist:
            site.energy = 0.
            site.gradient[:]  = 0.
        self.E = 0.
        self.comgrad[:] = 0.
        self.aagrad[:] = 0.

    def addSymmetryRotation(self, aa): #add rotational symmetry
        self.symmetrylist_rot.append( aa )
    #define functions to add new types of symmetries
    
    def getSymmetries(self, com, aa ):
        """
        a generator which iteratively returns the absolute xyz coordinates
        of the molecule subject to it's symmetries  
        
        com: the center of mass coords of the molecule
        
        aa: the orientation of the molecule in angle-axis  
        """
        com = np.array(com)
        xyz = np.zeros(3*self.nsites)
        rmat = rot.aa2mx(aa) #rotation matrix
        #first yield the unaltered molecule
        xyz = self.getxyz_rmat(rmat, com=com)
        yield xyz, aa
        #now loop through the symmetries
        for p in self.symmetrylist_rot:
            #combine the two rotations into one
            rmat_comb = np.dot( rmat, rot.aa2mx(p) )
            xyz = self.getxyz_rmat(rmat_comb, com=com)
            newaa = rot.rotate_aa(p, aa)
            #print rmat_comb
            #print rot.aa2mx( newaa)
            yield xyz, newaa
        #return other symmetries here
                
            
            

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
    
    otp.addSymmetryRotation( np.array([ 0., np.pi, 0.]))
    otp.correctCoM()
    return otp

def dumbbell(sig1 = 0.35, sig2 = 0.65):
    pos1 = [0., 0., 0.5]
    pos2 = [0., 0., -0.5]
    dbel = Molecule()
    dbel.insert_site(0, pos1)
    dbel.insert_site(1, pos2)
    dbel.correctCoM()
    
    from potentials.lj import LJ
    lj1 = LJ(sig = 2.*sig1)
    lj2 = LJ(sig = 2.*sig2)
    lj12 = LJ(sig = sig1+sig2)
    interaction_matrix = [ [lj1, lj12], [lj12, lj2] ]
    return dbel, interaction_matrix
    
    


def test_molecule():
    otp = setupLWOTP()

    xyz = otp.getxyz()
    from printing.print_atoms_xyz import printAtomsXYZ as printxyz
    import sys
    #with open("out.xyz", "w") as fout:
    printxyz(sys.stdout, xyz)
    
    aa = np.array([.2, .3, .4])
    for xyz, aanew in otp.getSymmetries( np.zeros(3), aa):
        printxyz(sys.stdout, xyz, line2="symmetry")
        xyz = otp.getxyz(aa = aanew)
        printxyz(sys.stdout, xyz, line2="symmetry from returned aa")




if __name__ == "__main__":
    test_molecule()
    
    aa1 = rot.random_aa()
    aa2 = rot.random_aa()
    rmat1 = rot.aa2mx(aa1)
    rmat2 = rot.aa2mx(aa2)
    
    rmat21 = np.dot(rmat2, rmat1)
    aa21 = rot.rotate_aa(aa1, aa2)
    rmat21aa = rot.aa2mx(aa21)
    print rmat21
    print rmat21aa
    print abs(rmat21 - rmat21aa) < 1e-12
