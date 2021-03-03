from __future__ import print_function
import string

__all__ = ["readAmberParam"]


class readAmberParam:
    """Extract info from coords.prmtop
    
    prmtop:
        string with the full coords.prmtop file 
    
    bondConn:
        list of 2-ples such that bondConn[i][0] and bondConn[i][1] are the 
        atom numbers of the i-th bond; atom numbers start at 1   
    
    TODO: deprecate -- use OpenMM parser instead (11 Jan 13) 
                
    """

    def __init__(self):
        prmtop = open("coords.prmtop").read()
        self.prmtopLines = string.split(prmtop, "\n")
        self.bondConn = []

    def populateBondConn(self):

        bondBlock = []  # list of numbers between q0 and q1
        tt = 0

        # # bonds with hydrogens
        #        # q0 and q1 are line numbers
        q0 = self.prmtopLines.index("%FLAG BONDS_INC_HYDROGEN                                                        ")
        q1 = self.prmtopLines.index(
            "%FLAG BONDS_WITHOUT_HYDROGEN                                                    ")  #

        for ind in range(q0 + 2, q1):
            for j in string.split(self.prmtopLines[ind]):
                bondBlock.append(int(j))
                tt += 1

        print('bonds containing hydrogen read')

        # bonds without hydrogen         
        q0 = self.prmtopLines.index("%FLAG BONDS_WITHOUT_HYDROGEN                                                    ")
        q1 = self.prmtopLines.index("%FLAG ANGLES_INC_HYDROGEN                                                       ")

        for ind in range(q0 + 2, q1):
            for j in string.split(self.prmtopLines[ind]):
                bondBlock.append(int(j))
                tt += 1

        print('bonds without hydrogen read')

        # total number of bonds 
        Nbonds = tt / 3

        # now populate bondConn 
        for i in range(Nbonds):
            ta1 = bondBlock[3 * i]
            ta2 = bondBlock[3 * i + 1]
            # good idea to sort it 
            if ta1 > ta2:
                self.bondConn.append((ta2 / 3 + 1, ta1 / 3 + 1))
            else:
                self.bondConn.append((ta1 / 3 + 1, ta2 / 3 + 1))

        print('bond connectivity read')

        # atom index in bondList starts from 1 

        # %FLAG BONDS_INC_HYDROGEN
        # %FORMAT(10I8)  (IBH(i),JBH(i),ICBH(i), i=1,NBONH)
        #  IBH    : atom involved in bond "i", bond contains hydrogen
        #  JBH    : atom involved in bond "i", bond contains hydrogen
        #  ICBH   : index into parameter arrays RK and REQ

        # the atom numbers in the following arrays that describe bonds, angles, and #dihedrals are coordinate array indexes for runtime speed. The true atom number #equals the absolute value of the number divided by three, plus one.

    def printBondConn(self):
        print('printing list of bonds')
        print(len(self.bondConn))

        for i in self.bondConn:
            print(i)


if __name__ == "__main__":
    mol = readAmberParam()
    mol.populateBondConn()
    mol.printBondConn()
    
    
    
    

