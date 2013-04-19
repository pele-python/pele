import numpy as np
import networkx as nx

mass_common_elements = {  1: "H",
                         12: "C",
                         14: "N",
                         16: "O",
                         32: "S" }

mass_all_elements = {}

class Atom(object):
    """
    A basic atom type which has a name, element, charge and LJ parameters of the form:
    
        V_lj(r) = A/r**12 - B/r**6
        
    Where V_lj(r) is the energy, A is lj_a and B is lj_b below.
    
    The atom types also maintain a simple list of their bonded atoms for easily 
    traversing their graph representation. 
    """
    def __init__(self, name, mass, charge, lj_a, lj_b):
        self.name = name
        self.mass = mass
        self.charge = charge
        self.set_element(mass)
        self.lj_a = lj_a
        self.lj_b = lj_b
        self.bonded = []
    def set_element(self, mass):
        """
        Defines the element, based on the mass. First, we look up in a dictionary of
        common (organic elements) masses. If that doesn't work, look up in the much
        longer list of all the elements.
        """
        self.int_mass = int(round(self.mass))
        try:
            self.element = mass_common_elements[self.int_mass]
        except KeyError:
            self.element = mass_all_elements[self.int_mass]
    def get_lj_sigma_epsilon(self):
        """
        This calculates sigma and epsilon, based on the LJ A and B coefficients
        described above and returns them as a tuple.
        
            s = (A/B)**(1/6)
            e = B**2/4A 
        """
        sigma   = np.reciprocal(np.power(self.lj_a/self.lj_b, 6.0))
        epsilon = self.lj_b * self.lj_b / (4.0 * self.lj_a)
        return sigma, epsilon 
    def add_bond(self, other_atom):
        """ Defines a bond between two Atoms. """
        self.bonded.append(other_atom)
        other_atom.bonded.append(self)
    def del_bond(self, other_atom):
        """ Deletes a bond between two Atoms. """
        if other_atom in self.bonded:
            self.bonded.remove(other_atom)
            other_atom.bonded.remove(self)
    def bond_count(self):
        """ Returns the number of bonds an Atom has (bond order is ignored). """
        return len(self.bonded)