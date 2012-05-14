import numpy as np
import molecule
import rotations as rot
import itertools
from potentials.potential import potential
import copy

class InteractionList:
    def __init__(self, interaction, typepair):
        self.interaction = interaction
        self.ilist = []
        self.types = [typepair]
        
    def addPair(self, ij ):
        self.ilist.append( ij )

    def addTypePair(self, i, j):
        self.types.append( (i,j) )

    def getNPilist(self):
        try:
            return self.npilist
        except:
            #convert ilist to np array
            self.npilist = np.array(self.ilist)
            return self.npilist
        
#def buildInteractionLists(typelist, interaction_matrix):
#    type2list = dict()
#    for i in typelist:
#        for j in typelist:

class RBSandbox(potential):
    """
    Defines a system of rigid body molecules
    """
    def __init__(self, molecule_list, interaction_matrix, use_interaction_lists = False):
        """
        molecule_list:  a list of Molecule objects that define the system

        interaction_matrix:  a matrix that defines the interactions between site types.
            inter = interaction_matrix[site1.type][site2.type] is the interaction between site1 and site2
            inter.getEnergy(coords) returns the energy of the interaction
            
        use_interaction_lists: indicates the potentials can accept a list of interacting particles
            This makes the calculation of energy, gradient much faster. 
        """
        self.molecule_list = [copy.deepcopy(mol) for mol in molecule_list]
        self.interaction_matrix = interaction_matrix
        self.nmol = len(self.molecule_list)

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

        self.buildInteractionLists()
        
        self.oldcoords = np.zeros(3*2*self.nmol)

    def getxyz(self, coords):
        """convert center of mass + angle-axis coords into xyz coordinates"""
        self.update_coords(coords)
        xyz = np.zeros(self.nsites*3)
        isite = 0
        for mol in self.molecule_list:
            for site in mol.sitelist:
                xyz[isite*3 : isite*3 + 3] = site.abs_position
                isite += 1
        return xyz

    def molmolEnergy(self, mol1, mol2):
        """
        return the energy of interaction between mol1 and mol2.

        don't update the molecules
        """
        E = 0.
        coords2 = np.zeros(6, np.float64)
        for site1 in mol1.sitelist:
            coords2[0:3] = site1.abs_position
            for site2 in mol2.sitelist:
                coords2[3:6] = site2.abs_position
                #print coords2
                interaction = self.interaction_matrix[site1.type][site2.type]
                dE = interaction.getEnergy(coords2)
                E += dE
        return E


    def getEnergy(self, coords):
        self.update_coords(coords)
        self.zeroEnergyGrad()
        E = 0.
        for imol1 in range(self.nmol):
            mol1 = self.molecule_list[imol1]
            for imol2 in range(0,imol1):
                mol2 = self.molecule_list[imol2]
                E += self.molmolEnergy(mol1, mol2)
        return E



    def coords_compare(self, coords):
        """ return true if coords is the same as oldcoords"""
        return all(coords == self.oldcoords)
        #maxdif = np.max( np.abs(coords - self.oldcoords))
        #return maxdif < 1e-16
            #print "maxdif", maxdif, coords[0:3]
            #return True
        #return False

    def update_coords(self, coords):
        """
        update the com and angle-axis coords and dependents on all molecules

        only do it if coords is different.
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
            mol.zeroEnergyGrad( )

    def getEnergyGradient(self, coords):
        self.update_coords(coords)
        self.zeroEnergyGrad()
        sitexyz = self.getxyz(coords)
        Etot, gradsite = self.siteEnergyGradient(sitexyz)
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
        return Etot, grad


    def updateSiteEnergyGradient(self, xyz, grad, ilist):
        pot = ilist.interaction
        try:
            e, g = pot.getEnergyGradientList(xyz, ilist.getNPilist())
            grad += g #this is a bad way of adding in the gradients. 
            return e
        except:
            print "using slow interaction lists potential"
            Etot = 0.
            coords2 = np.zeros(6, np.float64)
            for i,j in ilist.ilist:
                coords2[:3] = xyz[i*3:i*3+3]
                coords2[3:] = xyz[j*3:j*3+3]
                E, g = pot.getEnergyGradient(coords2)
                Etot += E
                grad[i*3:i*3+3] += g[:3]
                grad[j*3:j*3+3] += g[3:]
            return Etot

    def siteEnergyGradient(self, xyz):
        grad = np.zeros( self.nsites * 3 )
        Etot = 0.
        for ilist in self.ilists:
            Etot += self.updateSiteEnergyGradient( xyz, grad, ilist )
        return Etot, grad
    
    def buildInteractionLists(self):
        self.ilists = []
        type2ilist = dict()
        for mol1, mol2 in itertools.combinations( self.molecule_list, 2 ):
            for site1 in mol1.sitelist:
                t1 = site1.type
                i1 = site1.index
                for site2 in mol2.sitelist:
                    t2 = site2.type
                    i2 = site2.index
                    tp = (t1,t2)
                    tpi = (t2,t1)
                    if not tp in type2ilist:
                        # I should check here if some of the interactions are the same
                        inter = self.interaction_matrix[t1][t2]
                        type2ilist[tp] = InteractionList(inter, (t1,t2) )
                        self.ilists.append( type2ilist[tp] )
                        if t1 != t2:
                            type2ilist[tpi] = type2ilist[ tp ]
                    type2ilist[tp].addPair( (i1,i2) )
        



def test_symmetries(coords, mysys):
    nmol = mysys.nmol
    #the molecules now have their correct permutation.  
    #For each molecule, apply the symmetry operation which minimized the distance
    for i, mol in enumerate(mysys.molecule_list):
        comold = coords[3*i : 3*i + 3]
        kaa = 3*nmol + 3*i
        aaold = coords[kaa : kaa + 3]
        for xyz, aanew in mol.getSymmetries(comold, aaold):
            coords[kaa : kaa + 3] = aanew
            print "E symmetry", mysys.getEnergy(coords) 


def random_coords(nmol, nsites = None):
    if nsites == None: nsites = nmol
    coords = np.zeros(2*3*nmol)
    coords[0:nmol*3] = np.random.uniform(-1,1,[nmol*3]) * 1.3*(nsites)**(1./3)
    for i in range(nmol):
        k = nmol*3 + 3*i
        coords[k : k + 3] = rot.random_aa()
    return coords


def test_sandbox(nmol = 6):
    from numpy import sin, cos, pi
    import copy
    from potentials.lj import LJ

    #define the molecule types.
    #here use only one type, LWOTP
    otp = molecule.setupLWOTP()


    # define the interaction matrix for the system.
    # for LWOTP there is only one atom type, so this is trivial
    lj = LJ()
    interaction_matrix = [[lj]]


    #set up a list of molecules
    mols = [otp for i in range(nmol)]

    #set up the RBSandbox object
    mysys = RBSandbox(mols, interaction_matrix)
    nsites = mysys.nsites

    #get an initial set of coordinates
    comcoords = np.random.uniform(-1,1,[nmol*3]) * 1.3*(nsites)**(1./3)
    aacoords = np.array( [copy.copy(rot.random_aa()) for i in range(nmol)] )
    aacoords = aacoords.reshape(3*nmol)
    coords = np.zeros(2*3*nmol, np.float64)
    coords[0:3*nmol] = comcoords[:]
    coords[3*nmol:2*3*nmol] = aacoords[:]
    print "lencoords, len aacoords", len (coords), len(aacoords), len(comcoords)


    #save the initial set of coords
    printlist = []
    xyz = mysys.getxyz( coords )
    printlist.append( (xyz.copy(), "initial"))

    #calculate the initial energy
    Einit = mysys.getEnergy(coords)
    print "initial energy", Einit

    #test the gradient
    numericV = mysys.NumericalDerivative(coords, 1e-10)
    numericV = numericV.copy()
    print "numeric V", numericV

    E, V = mysys.getEnergyGradient(coords)
    print "energy from gradient", E, "difference", E - Einit
    #print "analytic gradient", V
    maxgrad_relative = np.max(np.abs(V-numericV)/np.abs(V))
    maxgraddiff = np.max(np.abs(V - numericV))
    print "max error in gradient", maxgraddiff, "max relative", maxgrad_relative

    #do a quench to make sure everything is working well
    import optimize.quench as quench
    coords, E, rms, funcalls = quench.quench(coords, mysys.getEnergyGradient, iprint=-1)
    print "postquench E", E, "rms", rms, "funtion calls", funcalls
    xyz = mysys.getxyz(coords )
    printlist.append( (xyz.copy(), "post quench"))


    #print the saved coords
    fname = "otp.xyz"
    print "saving xyz coords to", fname
    from printing.print_atoms_xyz import printAtomsXYZ as printxyz
    with open(fname, "w") as fout:
        for xyz, line2 in printlist:
            printxyz( fout, xyz, line2=line2, atom_type=["N", "O", "O"])
            
    
    test_symmetries(coords, mysys)

def test_sandbox_dumbbell():
    print "***********************************"
    print "testing dumbbell molecule"
    print "***********************************"

    dbel, imatrix = molecule.dumbbell()
    nmol = 5
    print imatrix[0][0] == imatrix[0][1]
    print imatrix[1][0] == imatrix[0][1]
    
    #set up a list of molecules
    mols = [dbel for i in range(nmol)]

    #set up the RBSandbox object
    mysys = RBSandbox(mols, imatrix)
    nsites = mysys.nsites

    coords = random_coords(nmol, nsites)
    
    #calculate the initial energy
    Einit = mysys.getEnergy(coords)
    print "initial energy", Einit

    #test the gradient
    numericV = mysys.NumericalDerivative(coords, 1e-10)
    numericV = numericV.copy()
    print "numeric V", numericV

    E, V = mysys.getEnergyGradient(coords)
    print "energy from gradient", E, "difference", E - Einit
    #print "analytic gradient", V
    maxgrad_relative = np.max(np.abs(V-numericV)/np.abs(V))
    maxgraddiff = np.max(np.abs(V - numericV))
    print "max error in gradient", maxgraddiff, "max relative", maxgrad_relative

    #do a quench to make sure everything is working well
    import optimize.quench as quench
    coords, E, rms, funcalls = quench.quench(coords, mysys.getEnergyGradient, iprint=-1)
    print "postquench E", E, "rms", rms, "funtion calls", funcalls

if __name__ == "__main__":
    test_sandbox()
    test_sandbox_dumbbell()

