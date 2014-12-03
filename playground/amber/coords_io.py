import numpy as np
import os

def write_inpcrd(coords, inpcrd_filename, name=None):
    """
    Writes coords to the file inpcrd_filename in the .inpcrd format.
    """
    # Open the output file.
    with open(inpcrd_filename, 'w') as inpcrd_file:
        print name
        try:
            # Print the name.
            inpcrd_file.write(name)
        except TypeError:
            # Name is not the appropriate type (e.g. None)
            pass
        # Print OS-specific newline character.
        inpcrd_file.write(os.linesep)
        # Print the number of atoms.
        inpcrd_file.write('{:5d}'.format(len(coords) / 3))
        inpcrd_file.write(os.linesep)
        # Now loop through coords.
        for index in range(0, len(coords), 6):
            coords_pair = coords[index:index+6]
            # For each coords_pair, there should be either three or six elements.
            # If there are three elements, the first should raise IndexError.
            try:
                inpcrd_file.write(''.join(['{:12.7f}']*6).format(*coords_pair))
            except IndexError:
                inpcrd_file.write(''.join(['{:12.7f}']*3).format(*coords_pair))
            inpcrd_file.write(os.linesep)

if __name__ == '__main__':
    # Test with odd number of atoms
    coords = np.random.random(15) - 0.5
    write_inpcrd(coords, 'my_test_file_15', '5_ATOMS')
    # Test with even number of atoms
    coords = np.random.random(18) - 0.5
    write_inpcrd(coords, 'my_test_file_18', '6_ATOMS')
