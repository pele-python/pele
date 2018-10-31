from pele.potentials import BasePotential
import numpy as np

__all__ = ["GMINPotential"]


class GMINPotential(BasePotential):  # pragma: no cover
    """
    Interface to fortran GMIN potential

    Potentials implemented in GMIN can be called from python if GMIN is compiled with the flag WITH_PYTHON enabled. This creates
    python modules (dynamic libraries). However, the interface is still very rough and GMINPotential provides a wrapper for
    easy access to GMIN.

    The imported GMIN module requires a data file to be present in the current directory. All parameters except for the ones
    responsible to setup the potential will be ignored and can be skipped. The first call after importing the module should be
    initialize.

    Attributes
    ----------

    GMIN :
        reference to the gmin module

    Examples
    --------

    The following example imports the GMIN python interface and evaluates the energy

    >>> import gmin_
    >>>
    >>> gmin_.initialize() # finish gmin initialization
    >>> pot = GMINPotential(gmin_)
    >>>
    >>> coords = pot.getCoords()
    >>> pot.getEnergy(coords)
    """

    def __init__(self, GMIN):
        """
        Constructor
        """
        self.GMIN = GMIN
        self.ncalls = 0

    def getEnergy(self, coords):
        self.ncalls += 1
        return self.GMIN.getEnergy(coords)

    def getEnergyGradient(self, coords):
        self.ncalls += 1
        grad = np.zeros(3 * self.GMIN.getNAtoms())
        E = self.GMIN.getEnergyGradient(coords, grad)
        return E, grad[0:coords.size]

    def getCoords(self):
        coords = np.zeros(self.GMIN.getDOF())
        self.GMIN.getCoords(coords)
        return coords

