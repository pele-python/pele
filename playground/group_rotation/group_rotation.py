from pele.takestep import TakestepInterface
import playground.group_rotation.transforms as transforms
import numpy as np

class GroupRotation(TakestepInterface):
    def __init__(self, group_rotation_params):
        self.parameters = group_rotation_params
    def takeStep(self, coords, **kwargs):
        step = {}
        for group in self.parameters:
            step[group] = {}
            # First check whether to actually take the step, if the uniform
            # random variable is greater than probability, go to the next group
            if np.random.uniform() > self.parameters[group]["selection_probability"]:
                step[group]["happened"] = 0
                step[group]["angle"] = 0.0
                step[group]["transform"] = np.identity(4)
                step[group]["pre_coords"] = coords.copy()
                step[group]["post_coords"] = coords.copy()
                continue
            # Extract the parameters from the dictionary
            bond_atom_1   = self.parameters[group]["bond_atom_1"]
            bond_atom_2   = self.parameters[group]["bond_atom_2"]
            max_angle_mag = self.parameters[group]["max_angle_magnitude"]
            group_atoms   = self.parameters[group]["group_atoms"]
            # Now define some more convenience values
            min_angle = -1.0 * np.pi * max_angle_mag
            max_angle = +1.0 * np.pi * max_angle_mag
            coords_1_indices = [3 * bond_atom_1 + k for k in range(3)]
            coords_2_indices = [3 * bond_atom_2 + k for k in range(3)]
            coords_1 = coords[coords_1_indices]
            coords_2 = coords[coords_2_indices]
            # Calculate an angle through which to rotate
            angle = np.random.uniform(low = min_angle, high = max_angle)
            # Store the information about what type of step will be taken
            step[group]["happened"] = 1
            step[group]["angle"] = angle
            step[group]["pre_coords"] = coords.copy()
            # Calculate the bond axis about which to rotate
            bond_axis = coords_2 - coords_1
            # Create the relevant AffineTransform object
            trans_to_origin = transforms.translation(-1.0 * coords_1)
            trans_back      = transforms.translation(+1.0 * coords_1)
            rotation        = transforms.proper_rotation(bond_axis,
                                                         angle,
                                                         affine=True)
            overall         = trans_back * rotation * trans_to_origin
            # Create a subset of the coords_array, N.B. using advanced indexing,
            # rather than normal indexing returns a copy of the indexed sections,
            # not a view. Thus, you can't edit the array in-place with advanced
            # indexing.
            coords_indices = [3 * atom + k for atom in group_atoms for k in range(3)]
            group_atom_coords = coords[coords_indices]
            # Apply the transform to the subset of coordinates
            group_atom_coords = overall(group_atom_coords)
            # Write the transformed subset back into the coordinates array
            for group_atom_coord_index, coord_index in enumerate(coords_indices):
                coords[coord_index] = group_atom_coords[group_atom_coord_index]
            # Set post_coords for the step taken
            step[group]["transform"] = overall
            step[group]["post_coords"] = coords.copy()
        return step

def averages(parameters, coords, n_trials, attribute):
    group_rotation = GroupRotation(parameters)
    attribute_average = {}
    attribute_variance = {}
    for group in parameters:
        attribute_average[group] = 0.0
        attribute_variance[group] = 0.0
    for i in range(n_trials):
        step_result = group_rotation.takeStep(coords)
        for group in parameters:
            attribute_average[group] = ((i) * attribute_average[group] +
                                       step_result[group][attribute]) / float(i + 1)
            diff = step_result[group][attribute] - attribute_average[group]
            attribute_variance[group] = ((i) * attribute_variance[group] +
                                        diff**2) / float(i + 1)
    return attribute_average, attribute_variance

def measure_dihedral(coords, atoms):
    indices = [3 * atom + k for atom in atoms for k in range(3)]
    coords_0 = coords[indices[0:3]]
    coords_1 = coords[indices[3:6]]
    coords_2 = coords[indices[6:9]]
    coords_3 = coords[indices[9:12]]
    b1 = coords_1 - coords_0
    b2 = coords_2 - coords_1
    b3 = coords_3 - coords_2
    b2_b3 = np.cross(b2, b3)
    b1_b2 = np.cross(b1, b2)
    angle = np.arctan2(np.linalg.norm(b2) * np.dot(b1, b2_b3), np.dot(b1_b2, b2_b3) )
    return angle
    
if __name__ == "__main__":
    import pele.amber.read_amber as amber
    import playground.group_rotation.amino_acids as amino
    
    topology_data = amber.read_topology("/home/khs26/coords.prmtop")
    parsed = amber.create_atoms_and_residues(topology_data)
    test_params = amber.group_rotation_dict(parsed, amino.def_parameters)
    test_coords = np.array(amber.read_amber_coords("/home/khs26/coords.inpcrd"))
    testGR = GroupRotation(test_params)
    pre_coords = test_coords.copy()
    result = testGR.takeStep(test_coords)
    print test_coords-pre_coords
