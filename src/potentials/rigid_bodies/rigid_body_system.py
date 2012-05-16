import numpy as np
from potentials.potential import potential as basepotential
#from potentials.rigid_body_potential import RigidBodyPotential
import copy


class RigidBodySystem(basepotential):
    """
    Defines a system of rigid body molecules
    """
    def __init__(self, molecule_list, potential = None):
        """
        molecule_list:  a list of Molecule objects that define the system
        
        potential:      the class which calculates site energies and gradients
                        It can be attached later via the function setPotential
        """
        self.molecule_list = [copy.deepcopy(mol) for mol in molecule_list]
        self.nmol = len(self.molecule_list)
        if potential != None:
            self.potential = potential

        self.nsites = 0
        self.typelist = []
        self.nsites_cum = np.zeros(self.nmol, np.int)
        sitenum = 0
        for i, mol in enumerate(self.molecule_list):
            self.nsites_cum[i] = self.nsites
            self.nsites += mol.nsites
            for site in mol.sitelist:
                self.typelist.append( site.type )
                site.index = sitenum
                sitenum += 1
        
        self.oldcoords = np.zeros(3*2*self.nmol)

    def setPotential(self, potential):
        """
        attach or replace the potential object.
        """
        self.potential = potential
        
    def transformToXYZ(self, coords):
        """
        convert center of mass + angle-axis coords into xyz coordinates of all the sites
        """
        self.update_coords(coords)
        xyz = np.zeros(self.nsites*3)
        isite = 0
        for mol in self.molecule_list:
            for site in mol.sitelist:
                xyz[isite*3 : isite*3 + 3] = site.abs_position
                isite += 1
        return xyz
    
    def updateGradientsOld(self, gradsite, coords):
        self.zeroEnergyGrad()
        grad = np.zeros([2*self.nmol*3])
        nmol = self.nmol
        for i, mol in enumerate(self.molecule_list):
            for site in mol.sitelist:
                    isite = site.index
                    site.gradient = gradsite[isite*3:isite*3+3]
                    mol.comgrad += site.gradient
                    #do angle axis part
                    mol.aagrad += np.dot( site.drdp, site.gradient)
            grad[         3*i :          3*i + 3] = mol.comgrad
            grad[3*nmol + 3*i : 3*nmol + 3*i + 3] = mol.aagrad
        return grad
    
    def updateGradients(self, sitegrad, coords):
        grad = np.zeros([2*self.nmol*3])
        isite = 0
        for i, mol in enumerate(self.molecule_list):
            iaa = self.nmol*3 + i*3
            comgrad, aagrad = mol.getGradients( coords[iaa:iaa+3], sitegrad[isite:isite + mol.nsites*3] )
            grad[i*3 : i*3 + 3] = comgrad
            grad[iaa : iaa + 3] = aagrad
            isite += mol.nsites*3
        return grad
        


    def getEnergy(self, coords):
        xyz = self.transformToXYZ(coords)
        return self.potential.getEnergy(xyz)

    def getEnergyGradient(self, coords):
        xyz = self.transformToXYZ(coords)
        E, xyzgrad = self.potential.getEnergyGradient(xyz)
        grad = self.updateGradients(xyzgrad, coords)
        return E, grad


    def getxyz(self, coords):
        return self.transformToXYZ(coords)

    def coords_compare(self, coords):
        """ return true if coords is the same as oldcoords"""
        return all(coords == self.oldcoords)

    def update_coords(self, coords):
        """
        update the com and angle-axis coords and dependents on all molecules

        only do it if coords has changed.
        """
        #using coords_compare makes quench almost always fail for some reason.
        if self.coords_compare( coords):
            return
        self.oldcoords[:] = coords[:]
        nmol= self.nmol
        for imol, mol in enumerate(self.molecule_list):
            com = coords[         imol*3 :          imol*3 + 3]
            aa  = coords[3*nmol + imol*3 : 3*nmol + imol*3 + 3]
            mol.update_coords( com, aa )

    def zeroEnergyGrad(self):
        for mol in self.molecule_list:
            mol.zeroEnergyGrad()
    

if __name__ == "__main__":
    pass
