"""
tools for reading from and writing to .xyz files

.. currentmodule:: pele.utils.xyz

.. autosummary::
    :toctree: generated/
    
    read_xyz
    write_xyz
"""

import numpy as np
from itertools import cycle
from collections import namedtuple

__all__ = ["read_xyz", "write_xyz"]


def read_xyz(fin):
    """ read a xyz file from file handle

    Parameters
    ----------
    fin : file handle
        file to read from

    Returns
    -------
    fin : open file
    xyz : namedtuple
        returns a named tuple with coords, title and list of atomtypes.

    See Also
    --------
    write_xyz

    """
    natoms = int(fin.readline())
    title = fin.readline()[:-1]
    coords = np.zeros([natoms, 3], dtype="float64")
    atomtypes = []
    for x in coords:
        line = fin.readline().split()
        atomtypes.append(line[0])
        x[:] = list(map(float, line[1:4]))

    return namedtuple("XYZFile", ["coords", "title", "atomtypes"]) \
        (coords, title, atomtypes)


def write_xyz(fout, coords, title="", atomtypes=("A",)):
    """ write a xyz file from file handle

    Writes coordinates in xyz format. It uses atomtypes as names. The list is
    cycled if it contains less entries than there are coordinates,

    One can also directly write xyz data which was generated with read_xyz.

    >>> xx = read_xyz("in.xyz")
    >>> write_xyz(open("out.xyz", "w"), *xx)

    Parameters
    ----------
    fout : an open file
    coords : np.array
        array of coordinates
    title : title section, optional
        title for xyz file
    atomtypes : iteratable
        list of atomtypes.

    See Also
    --------
    read_xyz

    """
    fout.write("%d\n%s\n" % (coords.size / 3, title))
    for x, atomtype in zip(coords.reshape(-1, 3), cycle(atomtypes)):
        fout.write("%s %.18g %.18g %.18g\n" % (atomtype, x[0], x[1], x[2]))
 
