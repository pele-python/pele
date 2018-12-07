import numpy as np

from pele.potentials import BasePotential
from pele.potentials.fortran import FinSin as fortran_fs

class FinSin(BasePotential):
    """The Finnis-Sinclair potential for transition metals.

    Parameters
    ----------
           d : float
               d is a cutoff for the pseudo-density 
           A : float
               A sets the energy scale for the attractive term
        beta : float
               beta adjusts the repulsive term at short distance
           c : float
               c is a cutoff for the pairwise repulsive term
    c0,c1,c2 : float
               Additional fitting parameters

    Notes
    -----
    The default parameter set is for Molybdenum, with the 
    values fitted to bulk properties and also adequate for 
    clusters.
    
    Here is a selection of representative parameter values:
    [ Phil. Mag. A 50, 45 (1984); ibid. 53, 161 (1986). ]

    Metal:   Vanadium    Chromium  Molybdenum    Tungsten
    
        d=   3.692767    3.915720    4.114825    4.400224
        A=   2.010637    1.453418    1.887117    1.896373
     beta=        0.0         1.8         0.0         0.0
        c=        3.8         2.9        3.25        3.25
      c_0= -0.8816318  29.1429813  43.4475218  47.1346499
      c_1=  1.4907756 -23.3975027 -31.9332978 -33.7665655
      c_2= -0.3976370   4.7578297   6.0804249   6.2541999
    
    """
    def __init__(self, d=4.114825, A=1.887117, beta=0.0, c=3.25,
                 c0=43.4475218, c1=-31.9332978, c2=6.0804249):
        self.d = d
        self.A = A
        self.beta = beta
        self.c = c
        self.c0 = c0
        self.c1 = c1
        self.c2 = c2
        
    def getEnergy(self, coords):
        e, g = self.getEnergyGradient(coords)
        return e

    def getEnergyGradient(self, coords):
        g, e = fortran_fs.finsin(coords, self.d, self.A, self.beta,
                                 self.c, self.c0, self.c1, self.c2)
        return e, g

#
# testing only below here
#
def test_molybdenum():
    pot = FinSin()
    natoms = 20
    x = np.random.uniform(-2,2,[3*natoms])
    pot.getEnergy(x)
    pot.test_potential(x)

if __name__ == "__main__":
    test_molybdenum()

