import numpy as np
import pygmin.utils.rotations as rot
from pygmin.potentials.fortran.rmdrvt import rmdrvt as rotMatDeriv
from scipy import weave
from scipy.weave import converters


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
    def __init__(self, sitetype, position):
        self.type = sitetype
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
        self.drmat = np.zeros([3,3,3]) #three 3x3 matrices for the derivatives of the rotation matrix
        self.comgrad = np.zeros(3)  #gradient of the center of mass
        self.aagrad = np.zeros(3)   #gradient of the angle-axis coords
        
        self.symmetrylist_rot = [] #list of rotational symmetries
        #TODO: implement other types of symmetry



    def insert_site(self, sitetype, position):
        """
        add a site to the definition of the molecule
        """
        self.sitelist.append( Site(sitetype, position) )
        self.nsites += 1
        self.sitexyz_molframe = np.array( [ site.position for site in self.sitelist ] )
        #print np.shape(self.sitexyz_molframe)
    
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
        self.sitexyz_molframe = np.array( [ site.position for site in self.sitelist ] )


    def getxyz_rmat(self, rmat, com = np.zeros(3)):
        """return the xyz positions of all sites in the molecule-frame"""
        xyz = np.transpose( np.dot( rmat, np.transpose(self.sitexyz_molframe) ) )
        xyz += com 
        xyz = np.reshape(xyz, self.nsites*3)
        return xyz


    def getxyz(self, com=np.array([0.,0.,0.]), aa=np.array([0.,0.,1e-6]) ):
        """return the xyz positions of all sites in the molecule-frame"""
        mx = rot.aa2mx(aa)
        return self.getxyz_rmat(mx, com)

    def update_rot_mat(self, aa, do_derivatives = True):
        self.rotation_mat, self.drmat[0,:,:], self.drmat[1,:,:], self.drmat[2,:,:] = rotMatDeriv(aa, do_derivatives)

    def getGradientsWeave(self, aa, sitegrad_in, recalculate_rot_mat = True):
        """
        convert site gradients to com and aa gradients
        """
        sitegrad = np.reshape(sitegrad_in, [self.nsites,3] )
        #calculate rotation matrix and derivatives
        if recalculate_rot_mat:
            self.update_rot_mat(aa, True)
        drmat = self.drmat
        x = self.sitexyz_molframe
        y = sitegrad
        aagrad = self.aagrad
        aagrad[:] = 0.  
        nsites= self.nsites    
        #self.aagrad = np.sum( np.sum( drmat * np.sum(x[:,np.newaxis,:]*y[:,:,np.newaxis] , axis=0), axis=2), axis=1 )
        #self.aagrad = np.sum( np.sum( drmat * np.dot( np.transpose(y), x) , axis=2), axis=1 )
        #for k in range(3):
            #self.aagrad[k] = np.sum( drmat[k,:,:] * np.dot( np.transpose(y), x) )
        """
        #The above looks very complicated, but it boils down to a three way sum.
        #The below is the simplest way to write it, but it's very slow
        self.aagrad = np.zeros(3)
        for k in range(3):
            for i in range(3):
                for j in range(3):
                    for isite in range(self.nsites):
                        self.aagrad[k] += drmat[k,i,j]*self.sitexyz_molframe[isite,j]*sitegrad[isite,i]
        """
        code = """
        for (int k=0; k<3; ++k){
            for (int i=0; i<3; ++i){
                for (int j=0; j<3; ++j){
                    for (int isite=0; isite<nsites; ++isite){
                        aagrad(k) += drmat(k,i,j)*x(isite,j)*y(isite,i);
                    }
                 }
            }
        }
        """
        weave.inline(code, ["aagrad", "drmat", "x", "y", "nsites"], type_converters=converters.blitz, verbose=2)
        self.comgrad = np.sum( sitegrad, axis=0 ) 
        self.aagrad = aagrad
        return self.comgrad, self.aagrad


    def getGradients(self, aa, sitegrad_in, recalculate_rot_mat = True):
        """
        convert site gradients to com and aa gradients
        """
        return self.getGradientsWeave(aa, sitegrad_in, recalculate_rot_mat)
        sitegrad = np.reshape(sitegrad_in, [self.nsites,3] )
        #calculate rotation matrix and derivatives
        if recalculate_rot_mat:
            self.update_rot_mat(aa, True)
        drmat = self.drmat
        x = self.sitexyz_molframe
        y = sitegrad      
        #self.aagrad = np.sum( np.sum( drmat * np.sum(x[:,np.newaxis,:]*y[:,:,np.newaxis] , axis=0), axis=2), axis=1 )
        self.aagrad = np.sum( np.sum( drmat * np.dot( np.transpose(y), x) , axis=2), axis=1 )
        #for k in range(3):
            #self.aagrad[k] = np.sum( drmat[k,:,:] * np.dot( np.transpose(y), x) )
        """
        #The above looks very complicated, but it boils down to a three way sum.
        #The below is the simplest way to write it, but it's very slow
        self.aagrad = np.zeros(3)
        for k in range(3):
            for i in range(3):
                for j in range(3):
                    for isite in range(self.nsites):
                        self.aagrad[k] += drmat[k,i,j]*self.sitexyz_molframe[isite,j]*sitegrad[isite,i]
        """
        self.comgrad = np.sum( sitegrad, axis=0 ) 
        return self.comgrad, self.aagrad

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
    """
    return a dumbell molecule.  i.e. two lj sites attached rigidly
    """
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
