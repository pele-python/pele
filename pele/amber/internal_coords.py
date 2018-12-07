from __future__ import print_function
import networkx as nx

import pele.amber.read_amber as ra


def cart_to_internal(molecule):
    """ Converts the coordinates of the molecule from cartesian to internal
    coordinates. """
    # First decide on a starting atom
    source_node = molecule.atoms.nodes()[0]
    # Argh, dictionaries aren't sorted...
    completed_bonds = {source_node: []}
    for edge in nx.dfs_edges(molecule.atoms, source_node):
        completed_bonds[edge[1]] = [edge[0]]
    # completed_bonds[2] += [x for x in completed_bonds[1][1:] if x not in completed_bonds[2]]
    # completed_bonds[3] += [x for x in completed_bonds[2][1:] if x not in completed_bonds[3]]

    for entry in completed_bonds:
        print(entry, completed_bonds[entry])
        # completed.append(completed[-1] + [search_tree.neighbors(source_node)[0]])
        #completed.append(completed[-1] + [[x for x in search_tree.neighbors(completed[-1][-1]) if x not in completed[-1]][0]])
        #completed.append(completed[-1] + [[x for x in search_tree.neighbors(completed[-1][-1]) if x not in completed[-1]][0]])

        #print nx.single_source_shortest_path(molecule.atoms, source_node)


if __name__ == "__main__":
    topology_data = ra.read_topology("/home/khs26/coords.prmtop")
    mol = ra.create_molecule(topology_data)
    mol.read_coords("/home/khs26/coords.inpcrd")
    cart_to_internal(mol)
