"""
Class containing sanity checks on proteins     
    
Parameters
----------
top   : open mm topology object created from prmtop file as 
        from .simtk.openmm.app import AmberPrmtopFile
        prmtop = AmberPrmtopFile('../../examples/amber/coords.prmtop')
        top = prmtop.topology  
"""
from __future__ import print_function
from __future__ import absolute_import

import numpy as np

__all__ = ["sanity_check"]


class sanity_check():
    def __init__(self, top):

        """
        top = OpenMM topology

        """
        print('in sanity check init')
        self.topology = top

        from . import measure

        self.measure = measure.Measure()

        self.populate_CAneighborList()
        self.populate_peptideAtomList()

    def populate_peptideAtomList(self):
        listofC = []
        listofO = []
        listofN = []
        listofH = []

        for i in self.topology.atoms():
            if i.name == 'C':
                listofC.append(i.index)

            if i.name == 'O':
                listofO.append(i.index)

            if i.name == 'N':
                listofN.append(i.index)

            if i.name == 'H':
                listofH.append(i.index)

                # print listofC
        # print listofO     
        # print listofN     
        # print listofH     

        # atom numbers of peptide bond 
        self.peptideBondAtoms = []

        for i in listofC:
            if listofO.__contains__(i + 1) and listofN.__contains__(i + 2) and listofH.__contains__(i + 3):
                self.peptideBondAtoms.append([i, i + 1, i + 2, i + 3])

        print('\nPeptide bond atom numbers (C,O,N,H, in order):  ')
        for i in self.peptideBondAtoms:
            print(i)

    def populate_CAneighborList(self):
        listofCA = []
        listofC = []
        listofN = []
        listofCB = []

        for i in self.topology.atoms():
            if i.name == 'CA':
                listofCA.append(i.index)

            if i.name == 'C':
                listofC.append(i.index)

            if i.name == 'N':
                listofN.append(i.index)

            if i.name == 'CB':
                listofCB.append(i.index)

                # print listofCA
        # print listofC     
        # print listofN     
        # print listofCB     

        # atom numbers of peptide bond 
        self.CAneighborList = []

        for i in listofCA:
            # find atoms bonded to CA 
            neighborlist = []
            for b in self.topology.bonds():
                if b[0] == i:
                    neighborlist.append(b[1])
                if b[1] == i:
                    neighborlist.append(b[0])

                    # print 'atoms bonded to CA ',i, ' = ', neighborlist
            nn = [i]
            # append C (=O) 
            for n in neighborlist:
                if listofC.__contains__(n):
                    nn.append(n)

                    # append CB
            for n in neighborlist:
                if listofCB.__contains__(n):
                    nn.append(n)

                    # append N
            for n in neighborlist:
                if listofN.__contains__(n):
                    nn.append(n)

            self.CAneighborList.append(nn)

            # atoms numbers start at 0
        print('\nCA neighbors atom numbers (CA,C(=O),CB, N, in order):  ')
        for i in self.CAneighborList:
            print(i)

    def check_cistrans(self, coords):
        """ 
        Sanity check on the isomer state of peptide bonds   
        
        Returns False if the check fails i.e. if any of the peptide bond is CIS         
        
        """

        isTrans = True

        for i in self.peptideBondAtoms:
            atNum = i[0]
            rC = np.array([coords[3 * atNum], coords[3 * atNum + 1], coords[3 * atNum + 2]])
            atNum = i[1]
            rO = np.array([coords[3 * atNum], coords[3 * atNum + 1], coords[3 * atNum + 2]])
            atNum = i[2]
            rN = np.array([coords[3 * atNum], coords[3 * atNum + 1], coords[3 * atNum + 2]])
            atNum = i[3]
            rH = np.array([coords[3 * atNum], coords[3 * atNum + 1], coords[3 * atNum + 2]])

            # compute O-C-N-H torsion angle 
            rad, deg = self.measure.torsion(rO, rC, rN, rH)

            # print 'peptide torsion (deg) ', i, ' = ', deg 
            # check cis 
            if deg < 90 or deg > 270:
                isTrans = False
            # print 'CIS peptide bond between atoms ', i, ' torsion (deg) = ', deg

        return isTrans

    def check_CAchirality(self, coords):
        """ 
        Sanity check on the CA to check if it is L of D    
        
        Returns False if the check fails i.e. if any D-amino acid is present          
        
        """

        # print 'in check CA chirality'

        isL = True

        for i in self.CAneighborList:
            atNum = i[0]
            rCA = np.array([coords[3 * atNum], coords[3 * atNum + 1], coords[3 * atNum + 2]])
            atNum = i[1]
            rC = np.array([coords[3 * atNum], coords[3 * atNum + 1], coords[3 * atNum + 2]])
            atNum = i[2]
            rCB = np.array([coords[3 * atNum], coords[3 * atNum + 1], coords[3 * atNum + 2]])
            atNum = i[3]
            rN = np.array([coords[3 * atNum], coords[3 * atNum + 1], coords[3 * atNum + 2]])

            # compute improper torsion angle between C-CA-CB and CA-CB-N 
            rad, deg = self.measure.torsion(rC, rCA, rCB, rN)

            # check cis 
            if deg < 180:
                # this condition was found by inspection of structures todo   
                isL = False
            # print 'chiral state of CA atom ', i[0], ' is D'
            # print 'CA improper torsion (deg) ', i, ' = ', deg

        return isL


if __name__ == "__main__":

    # create topology from prmtop file 
    from .simtk.openmm.app import AmberPrmtopFile

    prmtop = AmberPrmtopFile('../../examples/amber/coords.prmtop')

    scheck = sanity_check(prmtop.topology)

    # get coordinates from a pdb file   
    from .simtk.openmm.app import pdbfile as openmmpdbReader
    from .simtk.unit import angstrom as openmm_angstrom

    pdb = openmmpdbReader.PDBFile('../../examples/amber/coords.pdb')  # todo: coords.pdb is hardcoded 

    coords = pdb.getPositions() / openmm_angstrom
    coords = np.reshape(np.transpose(coords), 3 * len(coords), 1)

    # test 
    if scheck.check_CAchirality(coords):
        print('\nCA chirality test passed (all L)')
    else:
        print('\nCA chirality test passed (atleast one D-amino acid)')

    if scheck.check_cistrans(coords):
        print('\npeptide cis-trans test passed (all trans)')
    else:
        print('\npeptide cis-trans test failed (atleast one cis)')

    from . import coords2pdb

    coords2pdb.coords2pdb(coords, prmtop.topology, 'temp.pdb')
    
    
    

