import re
import exceptions as exc
  
 
data_casts = {"a": str,
              "A": str,
              "e": float,
              "E": float,
              "i": int,
              "I": int}
 
class AmberTopology():
    def __init__(self, *args, **kwargs):
        pass
    def import_from_file(self, filename):
        with open(filename, "r") as topology_file:
            parsed_data = self.parse_file(topology_file)
        return parsed_data
    def parse_file(self, file_object):
        current_flag = None
        current_format = None
        current_data_length = None
        current_data_type = None
        self.topology_data = {}
        for line in file_object:
            # Remove whitespace and newlines
            line = line.rstrip()
            if line.startswith("%VERSION"): continue
            elif line.startswith("%FLAG"):
                current_flag = line[6:]
                self.topology_data[current_flag] = []
            elif line.startswith("%FORMAT"):
                current_format      = line[8:-1]
                current_data_length = int(re.split("[aEI\.\)]", current_format)[1])
                current_data_type   = re.findall("[a-zA-Z]", current_format)[0]
            else:
                split_line = [line[i:i+current_data_length] for i in range(0, len(line), current_data_length)]
                formatted_line = map(data_casts[current_data_type], split_line)
                self.topology_data[current_flag] += formatted_line
        return self.topology_data
    def bonds(self, by_type):
        if not hasattr(self, "topology_data"):
            raise exc.ValueError("Topology file doesn't contain any data.")
        self.bonds = []
        topology_bonds_list = topology_data["BONDS_INC_HYDROGEN"] + topology_data["BONDS_WITHOUT_HYDROGEN"]
        if by_type is True:
            while topology_bonds_list:
                atom_1    = topology_bonds_list.pop(0) / 3
                atom_2    = topology_bonds_list.pop(0) / 3
                bond_type = topology_bonds_list.pop(0)
                self.bonds.append((atom_1, atom_2, bond_type))
        else:
            while topology_bonds_list:
                atom_1    = topology_bonds_list.pop(0) / 3
                atom_2    = topology_bonds_list.pop(0) / 3
                bond_type = topology_bonds_list.pop(0)
                self.bonds.append((atom_1, atom_2))
        return self.bonds
    def angles(self, by_type):
        if not hasattr(self, "topology_data"):
            raise exc.ValueError("Topology file doesn't contain any data.")
        self.angles = []
        topology_angles_list = topology_data["ANGLES_INC_HYDROGEN"] + topology_data["ANGLES_WITHOUT_HYDROGEN"]
        if by_type is True:
            while topology_angles_list:
                atom_1    = topology_angles_list.pop(0) / 3
                atom_2    = topology_angles_list.pop(0) / 3
                atom_3    = topology_angles_list.pop(0) / 3
                bond_type = topology_angles_list.pop(0)
                self.angles.append((atom_1, atom_2, atom_3, bond_type))
        else:
            while topology_angles_list:
                atom_1    = topology_angles_list.pop(0) / 3
                atom_2    = topology_angles_list.pop(0) / 3
                atom_3    = topology_angles_list.pop(0) / 3
                bond_type = topology_angles_list.pop(0)
                self.angles.append((atom_1, atom_2, atom_3))
        return self.angles
    def dihedrals(self, by_type):
        if not hasattr(self, "topology_data"):
            raise exc.ValueError("Topology file doesn't contain any data.")
        self.dihedrals = []
        topology_dihedrals_list = topology_data["DIHEDRALS_INC_HYDROGEN"] + topology_data["DIHEDRALS_WITHOUT_HYDROGEN"]
        if by_type is True:
            while topology_dihedrals_list:
                atom_1    = topology_dihedrals_list.pop(0) / 3
                atom_2    = topology_dihedrals_list.pop(0) / 3
                atom_3    = topology_dihedrals_list.pop(0) / 3
                atom_4    = topology_dihedrals_list.pop(0) / 3
                bond_type = topology_dihedrals_list.pop(0)
                self.dihedrals.append((atom_1, atom_2, atom_3, atom_4, bond_type))
        else:
            while topology_dihedrals_list:
                atom_1    = topology_dihedrals_list.pop(0) / 3
                atom_2    = topology_dihedrals_list.pop(0) / 3
                atom_3    = topology_dihedrals_list.pop(0) / 3
                atom_4    = topology_dihedrals_list.pop(0) / 3
                bond_type = topology_dihedrals_list.pop(0)
                self.dihedrals.append((atom_1, atom_2, atom_3, atom_4))
        return self.dihedrals

#################################################### MAIN ####################################################

if __name__ == "__main__":
    amber_top = AmberTopology()
    topology_data = amber_top.import_from_file("/home/khs26/coords.prmtop")
#     for key in topology_data:
#         print key, topology_data[key]
    for bond in amber_top.bonds(by_type = True):
        print bond
    for angle in amber_top.angles(by_type = True):
        print angle
    for dihedral in amber_top.dihedrals(by_type = True):
        print dihedral
