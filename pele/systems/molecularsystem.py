from pele.systems import BaseSystem


class MolecularSystem(BaseSystem):
    """
    Representation for a molecular system, this system stores info about atoms, bonds,
    angles and torsions.
    
    It is possible to represent the molecule using a graph. However, this class is used
    to quickly and efficiently:
        - add/remove atoms, bonds, angles and torsions;
        - read/write PDB and other common formats;
        - interface between different formats of input files for CHARMM, AMBER etc.;
        - visualise molecular structures;
        - measure distances between structures.
    """


class Atom(object):
    """
    Representation of an Atom, object. Can have 
    """
