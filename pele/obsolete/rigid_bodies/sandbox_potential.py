import numpy as np
import itertools
from pele.potentials.potential import potential as basepotential
#from potentials.rigid_body_potential import RigidBodyPotential
        
        
class InteractionList:
    """
    a class to define and keep track of interaction lists
    """
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
        except AttributeError: #if self.npilist is not defined yet
            #convert ilist to np array, but do it only once.  It's quite slow
            self.npilist = np.array(self.ilist)
            return self.npilist


class RBSandboxPotential(basepotential):
    """
    This class builds a potential object for rigid bodies
    """
    def __init__(self, RB, interaction_matrix):
        """
        RB:  a RigidBodySystem object

        interaction_matrix:  a matrix that defines the interactions between site types.
            inter = interaction_matrix[site1.type][site2.type] is the potential class 
                    between site type site1 and site type site2
            inter.getEnergy(coords) returns the energy of the interaction
        """

        self.buildInteractionLists(RB, interaction_matrix)
        #for ilist in self.ilists:
        #    print len(ilist.ilist)


    def getEnergy(self, xyz):
        Etot = 0.
        for ilist in self.ilists:
            pot = ilist.interaction
            Etot += pot.getEnergyList(xyz, ilist.getNPilist() )
        return Etot


    def getEnergyGradient(self, xyz):
        """
        return the energy and site-gradient from all site-site interactions
        """
        grad = np.zeros( len(xyz) )
        Etot = 0.
        for ilist in self.ilists:
            #get contribution from interaction list ilist
            potential = ilist.interaction
            de, dg = potential.getEnergyGradientList( xyz, ilist.getNPilist() )
            Etot += de
            grad += dg #there may be more efficient ways of adding in the gradient contributions
        return Etot, grad
    
    def buildInteractionLists(self, RB, interaction_matrix ):
        """
        build interaction lists from site types and interaction_matrix
        
        This builds self.ilist, which is a list interaction lists.  
            An interaction list is a list of site pairs which interact with 
                the same potential
            There will be one interaction list for each potential type.
            e.g. i1, i2 = self.ilist[0][0] where i1 and i2 are the indices of 
                sites in the xyz representation 
            
        """
        self.ilists = []
        type2ilist = dict()
        for mol1, mol2 in itertools.combinations( RB.molecule_list, 2 ):
            for site1 in mol1.sitelist:
                t1 = site1.type
                i1 = site1.index
                for site2 in mol2.sitelist:
                    #find out what is the potential between site1 and site2
                    t2 = site2.type
                    i2 = site2.index
                    tp = (t1,t2)
                    tpi = (t2,t1)
                    if not tp in type2ilist:
                        # I should check here if some of the interactions are the same
                        inter = interaction_matrix[t1][t2]
                        type2ilist[tp] = InteractionList(inter, (t1,t2) )
                        self.ilists.append( type2ilist[tp] )
                        if t1 != t2: #make type2ilist symmetric
                            type2ilist[tpi] = type2ilist[ tp ]
                    type2ilist[tp].addPair( (i1,i2) )
        

if __name__ == "__main__":
    pass
