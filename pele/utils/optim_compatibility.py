"""
tools for reading and writing OPTIM input and output files
"""

import numpy as np

_id_count = 0

class UnboundMinimum(object):
    """
    a class to duplicate some of the functionality of the Minimum class
    """
    energy = None
    coords = None
    fvib = None
    pgorder = None
#    _id_count = 0
    def __init__(self, energy, coords):
        self.energy = energy
        self.coords = coords
        global _id_count
        self._id = _id_count
        _id_count += 1

    def __eq__(self, m):
        """m can be integer or Minima object"""
        assert self._id is not None
        if isinstance(m, UnboundMinimum):
            assert m._id is not None
            return self._id == m._id
        else:
            return self._id == m
        
    def __hash__(self):
        assert self._id is not None
        return self._id

def read_min_data(fname="min.data"):
    """
    return a list of minima with data read in from a min.data file
    
    Parameters
    ----------
    fname : str
        the file name of the min.data file
    
    Returns
    -------
    list of UnboundMinimum objects.
    """
    minima = []
    
    with open(fname, "r") as fin:
        for line in fin:
            sline = line.split()
            energy = float(sline[0])
            coords = np.array([0.])
            m = UnboundMinimum(energy, coords)
            
            m.fvib = float(sline[1])
            m.pgorder = int(sline[2])
            
            minima.append(m)
    return minima