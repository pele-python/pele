"""
Function to convert xyz coordinates to pdb 
"""
from __future__ import print_function

from .simtk.openmm.app import pdbfile as openmm_pdb
from .simtk.openmm import unit as openmm_unit
from .simtk.unit import angstrom as openmm_angstrom

__all__ = ["coords2pdb"]


def coords2pdb(coords, top, pdbfname):
    """
    coords   = 3*N x 1 numpy array 
    top      = open mm topology 
    pdbfname = output filename
    """

    fpointer = open(pdbfname, 'w')
    print("REMARK ", file=fpointer)

    coordinates = openmm_unit.Quantity(coords.reshape(top._numAtoms, 3), openmm_angstrom)

    openmm_pdb.PDBFile.writeModel(top, coordinates, fpointer)

    print("END", file=fpointer)

    fpointer.close()

    print('temp.pdb created')


#
# ct = 0
# for coords in coordslist:
#            ct = ct + 1 
#            coords = CoMToOrigin(coords.copy())
#            self.potential.copyToLocalCoords(coords) 
##            openmmpdb.PDBFile.writeFile(self.potential.prmtop.topology , self.potential.localCoords * openmm_angstrom , file=sys.stdout, modelIndex=1)
#            openmmpdb.PDBFile.writeModel(self.potential.prmtop.topology , self.potential.localCoords * openmm_angstrom , file=f, modelIndex=ct)
#                        
#        print "closing file"
#        f.flush()
#                
#        # load the molecule from the temporary file
#        pymol.cmd.load(fname)

