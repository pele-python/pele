#!/usr/bin/python

#####################################################
#
# TO DO
#
# Write group rotation script
# Work out bonding within the amino acid residue
# Implement BFS
# Consider making the standard rotations "special"
# Cycle detection?
#
#####################################################

import playground.group_rotation.search as search
import itertools
import re

amino_acids = [ "ARG",
                "HIE",
                "HID",
                "HIP",
                "LYS",
                "ASP",
                "GLU",
                "SER",
                "THR",
                "ASN",
                "GLN",
                "CYS",
                "GLY",
                "PRO",
                "ALA",
                "ILE",
                "LEU",
                "MET",
                "PHE",
                "TRP",
                "TYR",
                "VAL" ]

def_parameters = {}
# C, alpha-C bonds
for amino_acid in amino_acids:
    def_parameters[(amino_acid, ("C", "CA"))] = (0.1, 0.5)
# alpha-C, beta-C bonds
def_parameters[("ARG", ("CA", "CB"))] = (0.2, 0.5)
def_parameters[("HIE", ("CA", "CB"))] = (0.3, 0.5)
def_parameters[("HID", ("CA", "CB"))] = (0.3, 0.5)
def_parameters[("HIP", ("CA", "CB"))] = (0.3, 0.5)
def_parameters[("LYS", ("CA", "CB"))] = (0.2, 0.5)
def_parameters[("ASP", ("CA", "CB"))] = (0.5, 0.5)
def_parameters[("GLU", ("CA", "CB"))] = (0.3, 0.5)
def_parameters[("SER", ("CA", "CB"))] = (1.0, 0.5)
def_parameters[("THR", ("CA", "CB"))] = (1.0, 0.5)
def_parameters[("ASN", ("CA", "CB"))] = (0.5, 0.5)
def_parameters[("GLN", ("CA", "CB"))] = (0.3, 0.5)
def_parameters[("CYS", ("CA", "CB"))] = (1.0, 0.5)
def_parameters[("ALA", ("CA", "CB"))] = (1.0, 0.5)
def_parameters[("VAL", ("CA", "CB"))] = (1.0, 0.5)
def_parameters[("ILE", ("CA", "CB"))] = (0.5, 0.5)
def_parameters[("LEU", ("CA", "CB"))] = (0.5, 0.5)
def_parameters[("MET", ("CA", "CB"))] = (0.5, 0.5)
def_parameters[("PHE", ("CA", "CB"))] = (0.3, 0.5)
def_parameters[("TYR", ("CA", "CB"))] = (0.3, 0.5)
def_parameters[("TRP", ("CA", "CB"))] = (0.3, 0.5)
# beta-C, gamma-C bonds
def_parameters[("ARG", ("CB", "CG"))] = (0.3, 0.5)
def_parameters[("HIE", ("CB", "CG"))] = (0.5, 0.5)
def_parameters[("HID", ("CB", "CG"))] = (0.5, 0.5)
def_parameters[("HIP", ("CB", "CG"))] = (0.5, 0.5)
def_parameters[("LYS", ("CB", "CG"))] = (0.3, 0.5)
def_parameters[("ASP", ("CB", "CG"))] = (1.0, 0.5)
def_parameters[("GLU", ("CB", "CG"))] = (0.5, 0.5)
def_parameters[("ASN", ("CB", "CG"))] = (1.0, 0.5)
def_parameters[("GLN", ("CB", "CG"))] = (0.5, 0.5)
def_parameters[("ILE", ("CB", "CG1"))] = (1.0, 0.5)
def_parameters[("LEU", ("CB", "CG"))] = (1.0, 0.5)
def_parameters[("MET", ("CB", "CG"))] = (0.7, 0.5)
def_parameters[("PHE", ("CB", "CG"))] = (0.5, 0.5)
def_parameters[("TYR", ("CB", "CG"))] = (0.5, 0.5)
def_parameters[("TRP", ("CB", "CG"))] = (0.4, 0.5)
# gamma-C, delta-C bonds
def_parameters[("ARG", ("CG", "CD"))] = (0.5, 0.5)
def_parameters[("LYS", ("CG", "CD"))] = (0.5, 0.5)
def_parameters[("GLU", ("CG", "CD"))] = (1.0, 0.5)
def_parameters[("GLN", ("CG", "CD"))] = (1.0, 0.5)
def_parameters[("MET", ("CG", "SD"))] = (0.5, 0.5)
# delta, epsilon bonds
def_parameters[("ARG", ("CD", "NE"))] = (0.5, 0.5)
def_parameters[("LYS", ("CD", "CE"))] = (0.5, 0.5)

class Atom(object):
    """ Atom defined from the AMBER topology file. """
    def __init__(self, index, name, mass, amber_atom_type, charge):
        self.index = index
        self.name = name
        self.mass = mass
        self.amber_atom_type = amber_atom_type
        self.charge = charge
        self.residue = None
        self.set_element(mass)
        self.bonded = []
    def set_element(self, mass):
        """ Defines the element, based on the mass. """
        if 0.9 < mass and 1.1 > mass:
            self.element = "H"
        elif 11.9 < mass and 12.1 > mass:
            self.element = "C"
        elif 13.9 < mass and 14.1 > mass:
            self.element = "N"
        elif 15.9 < mass and 16.1 > mass:
            self.element = "O"
        elif 31.9 < mass and 32.1 > mass:
            self.element = "S"
    def set_residue(self, residue):
        """ Defines the Residue to which the Atom belongs. """
        self.residue = residue
    def add_bond(self, other_atom):
        """ Defines a bond between two Atoms. """
        self.bonded.append(other_atom)
        other_atom.bonded.append(self)
    def del_bond(self, other_atom):
        """ Deletes a bond between two Atoms. """
        if other_atom in bonded:
            self.bonded.remove(other_atom)
            other_atom.bonded.remove(self)
    def bond_count(self):
        """ Returns the number of bonds an Atom has (bond order is ignored). """
        return len(self.bonded)

class Residue(object):
    """ Residue defined from the AMBER topology file. """
    def __init__(self, index, name):
        self.index = index
        self.name = name
        self.atoms = []
    def add_atoms(self, atoms):
        """ 
        Adds Atoms in atoms to the Residue and sets the residue for those
        Atoms accordingly.
        """
        for atom in atoms:
            atom.set_residue(self)
            self.atoms.append(atom)

def split_len(seq, length):
    """ Returns a list containing the elements of seq, split into length-long chunks. """
    return [seq[i:i+length] for i in range(0, len(seq), length)]

def data_cast(datum, type):
    """
    Casts the datum as defined by type.

    a = string
    i = integer
    e = float (in exponent notation)
    """
    if type == "a" or type == "A":
        return str(datum)
    if type == "i" or type == "I":
        return int(datum)
    if type == "e" or type == "E":
        return float(datum)

def read_topology(filename):
    """ Reads a topology file, filename, and returns data blocks containing the topology data. """
    # Read the topology file to a series of lines. 
    with open(filename, "r") as input:
        prmtop_input_lines = []
        for line in input:
            line = line.strip("\n")
            if line.startswith('%VERSION'):
                continue
            if line.startswith('%FLAG'):
                prmtop_input_lines.append([line])
            else:
                prmtop_input_lines[-1].append(line)
    # Parse the input into data blocks, a list containing dictionaries with the following entries.
    # "flag": flag describing the data block (e.g. atom number, residue type);
    # "format": describes the format of the data block (e.g. 8a10 is 8 data of type char, length 10);
    # "data": contains the data, split into blocks of appropriate length and formatted to the correct type.
    prmtop_data_blocks = {}
    for input_line in prmtop_input_lines:
        # Get rid of %FLAG from front and strip whitespace.
        flag = input_line[0][6:].rstrip(" ")
        # Get rid of %FORMAT from front and strip ) and whitespace.
        format = input_line[1][8:].rstrip(") ")
        # Find length and data type.
        data_length = int(re.split("[aEI\.\)]",input_line[1])[1])
        data_type = re.findall("([a-zA-Z])", input_line[1][8:])[0]
        # Split the string into chunks of the appropriate size and cast to appropriate data type.
        split_string = split_len("".join(input_line[2:]), data_length)
        data = map(data_cast, split_string, itertools.repeat(data_type, len(split_string)))
        prmtop_data_blocks[flag] = data
    return prmtop_data_blocks

def create_atoms_and_residues(topology_data):
    """
    This function takes a topology data set (topology_data) and returns a set of Atoms and Residues
    using the types defined above.
    
    It also defines bonds between any bonded Atoms according to the topology file.
    """
    # First create lists of residues and atoms from the topology data. 
    atoms = map(Atom,
                range(1, (topology_data["POINTERS"][0] + 1)),      # The first element of POINTERS is the atom count.
                topology_data["ATOM_NAME"],
                topology_data["MASS"],
                topology_data["AMBER_ATOM_TYPE"],
                topology_data["CHARGE"])
    residues = map(Residue,
                   range(1, (topology_data["POINTERS"][11] + 1)),  # The 11th element of POINTERS is the residue count.
                   topology_data["RESIDUE_LABEL"])
    # Next put the appropriate atoms into the appropriate residues.  AMBER specifies the first atom index of
    # each residue (i.e. no end, hence the exception).
    residue_indices = topology_data["RESIDUE_POINTER"]
    for i, residue in enumerate(residues):
        try:
            start = residue_indices[i] - 1
            end = residue_indices[i+1] - 1
        except IndexError, ex:            
            end = None
        residue.add_atoms(atoms[start:end])
    # Go through the BONDS_INC_HYDROGEN and BONDS_WITHOUT_HYDROGEN lists to extract lists of bonded atoms.
    # Index is i / 3 + 1, because AMBER (in its infinite wisdom) still uses coordinate indices for runtime speed.
    first_atoms  = map(lambda x: (x / 3) + 1,
                       topology_data["BONDS_INC_HYDROGEN"][0::3] + topology_data["BONDS_WITHOUT_HYDROGEN"][0::3])
    second_atoms = map(lambda x: (x / 3) + 1,
                       topology_data["BONDS_INC_HYDROGEN"][1::3] + topology_data["BONDS_WITHOUT_HYDROGEN"][1::3])
    bond_list = zip(first_atoms, second_atoms)
    # Now go through bond list and create bonds between the relevant atoms.
    for bond in bond_list:
        atoms[(bond[0] - 1)].add_bond(atoms[(bond[1] - 1)])
    return (atoms, residues)

def group_rotation(atoms, residues, params, filename):
    """
    Defines group rotation for bonds between atoms defined in params.  This function DOESN'T check
    whether the atoms are actually bonded at any point (so make sure you're referencing appropriate
    atom pairs.

    atoms and residues are the atoms and residues lists returned from create_atoms_and_residues()

    params has the format:

    params[(residue, (atom_x, atom_y))] = (rotation_scale, probability)
    e.g. params[("ARG", ("CA", "CB"))] = (0.5, 0.2) rotates about the ARG CA-CB bond at most 90 degs each 
         time, with probability 0.2.
    """
    # Loop through the parameters, creating GROUPROTATION information for any bonds described by an individual
    # parameter.
    with open(filename, "w") as output:
        for param in params.keys():
            for residue in residues:
                if residue.name == param[0].ljust(4):
                    rotated_atoms = get_rotated_atoms(atoms, residues, (residue, param[1]))
                    # Create an identifiable group name consisting of [index]_[res_name]_[atom_1]_[atom_2]
                    group_name = "_".join((str(residue.index), 
                                           residue.name.strip(), 
                                           rotated_atoms[0].name.strip(),
                                           rotated_atoms[1].name.strip()))
                    output.write(" ".join(("GROUP", 
                                           group_name,
                                           str(rotated_atoms[0].index),
                                           str(rotated_atoms[1].index),
                                           str(len(rotated_atoms[2])),
                                           str(params[param][0]),
                                           str(params[param][1]),
                                           "\n")))
                    for atom in rotated_atoms[2]:
                        output.write(str(atom.index) + "\n")
                    output.write("\n")

def get_rotated_atoms(atoms, residues, bond):
    """
    Returns a list of the atoms in the rotating group for the GROUPROTATION script.
    
    bond has the format:

    bond = (residue, (atom_x, atom_y))
    """
    blocked_atoms = []
    for res_atom in bond[0].atoms:
        # Set the blocked atom (i.e. search for everything on one side of the bond) to
        # bond[1][0].
        if res_atom.name == bond[1][0].ljust(4):
            blocked_atoms.append(res_atom)
        # Set the root atom for the search to bond[1][1]
        if res_atom.name == bond[1][1].ljust(4):
            root_atom = res_atom
    return (blocked_atoms[0], root_atom, search.find_children(root_atom, blocked_atoms))

#for atom in atoms:
#    print atom.name, atom.index, (atom.residue).name, atom.bond_count()
#    if atom.name in ("C   ",):
#        blocked.append(atom)
#    if atom.name in ("CA  ",) and atom.residue.index == 2:
#        root = atom
#
#print map(lambda x: (x.name, x.index), search.find_children(root, blocked))
